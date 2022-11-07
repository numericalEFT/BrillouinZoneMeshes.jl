using BrillouinZoneMeshes
using BrillouinZoneMeshes.BZMeshes
using BrillouinZoneMeshes.BZMeshes.Coordinates
using BrillouinZoneMeshes.CompositeGrids
using BrillouinZoneMeshes.BaseMesh
using BrillouinZoneMeshes.AbstractMeshes
using BrillouinZoneMeshes.LinearAlgebra
using BrillouinZoneMeshes.StaticArrays

using Roots
using Test

struct RescaledGrid{T,GT,FT,IFT} <: AbstractGrid{T}
    # rescaled grid with function func()
    bound::SVector{2,T}
    size::Int
    grid::Vector{T} # grid == func.(basegrid.grid)

    #additional info
    basegrid::GT
    func::FT
    invfunc::IFT
    function RescaledGrid(basegrid::GT,
        func::FT, invfunc::IFT) where {GT,FT,IFT}

        T = eltype(GT)
        bound = func.(basegrid.bound)
        size = length(basegrid)
        grid = func.(basegrid.grid)

        return new{T,GT,FT,IFT}(bound, size, grid, basegrid, func, invfunc)
    end
end

# many interface could simply inherit from AbstractGrid
# the following need new implementation
Base.floor(grid::RescaledGrid, x) = floor(grid.basegrid, grid.invfunc(x))
CompositeGrids.Interp.locate(grid::RescaledGrid, x) = CompositeGrids.Interp.locate(grid.basegrid, grid.invfunc(x))

function radial_rescale(; grid::AbstractGrid, DIM::Int)
    if DIM == 2
        func = sqrt
        invfunc = x -> x^2
    elseif DIM == 3
        func = cbrt
        invfunc = x -> x^3
    else
        error("DIM=$DIM not implemented!")
    end

    return RescaledGrid(grid, func, invfunc)
end

function find_kFermi(dispersion, angle; kinit=0.0)
    if length(angle) == 1
        return find_zero(k -> dispersion(BZMeshes._polar2cart(Polar(k, angle...))), kinit)
    elseif length(angle) == 2
        return find_zero(k -> dispersion(BZMeshes._spherical2cart(Spherical(k, angle...))), kinit)
    else
        error("dimension $(length(angle)+1) not implemented!")
    end
end

function RescaledLogDensedGrid(type, bound, densepoints, Nlog, minterval, Nbase, DIM)
    if DIM == 2
        func = sqrt
        invfunc = x -> x^2
    elseif DIM == 3
        func = cbrt
        invfunc = x -> x^3
    else
        error("DIM=$DIM not implemented!")
    end
    rbound = invfunc.(bound)
    rdp = invfunc.(densepoints)
    g = CompositeGrid.LogDensedGrid(type, rbound, rdp, Nlog, minterval, Nbase)
    return radial_rescale(grid=g, DIM=DIM)
end

function kF_densed_kgrids(; dispersion,
    anglemesh,
    bound,
    basegridtype=:cheb,
    Nloggrid=3,
    minterval=0.01,
    Nbasegrid=2,
    DIM=2)
    # assume dispersion==0 has one root for each angle
    k1 = find_kFermi(dispersion, anglemesh[1])
    # g1 = CompositeGrid.LogDensedGrid(basegridtype, bound, [k1,], Nloggrid, minterval, Nbasegrid)
    g1 = RescaledLogDensedGrid(basegridtype, bound, [k1,], Nloggrid, minterval, Nbasegrid, DIM)
    grids = [g1,]
    for i in 2:length(anglemesh)
        kF = find_kFermi(dispersion, anglemesh[i])
        # g = CompositeGrid.LogDensedGrid(basegridtype, bound, [kF,], Nloggrid, minterval, Nbasegrid)
        g = RescaledLogDensedGrid(basegridtype, bound, [kF,], Nloggrid, minterval, Nbasegrid, DIM)
        push!(grids, g)
    end
    return grids
end


function BZMeshes.PolarMesh(; dispersion, anglemesh, br, kmax,
    kwargs...)

    DIM = size(br.lattice, 1)
    bound = [0.0, kmax]
    grids = kF_densed_kgrids(dispersion=dispersion, anglemesh=anglemesh, bound=bound, DIM=DIM, kwargs...)
    cm = CompositeMesh(anglemesh, grids)
    pm = PolarMesh(br, cm)
    return pm
end

@testset "PolarMesh Generator" begin
    @testset "2D" begin
        # given dispersion function accept a k in cartesian
        # goal is to find k_F at direction specified by angle
        dispersion(k) = dot(k, k) - 1.0

        # 2d
        N = 10
        bound = [-π, π]
        theta = SimpleGrid.Uniform(bound, N; isperiodic=true)

        k_F_previous = 0.0
        for θ in theta
            f(k) = dispersion(BZMeshes._polar2cart(Polar(k, θ)))
            k_F = find_zero(f, k_F_previous)
            @test k_F ≈ 1.0
            @test find_kFermi(dispersion, θ; kinit=k_F_previous) ≈ 1.0
            k_F_previous = k_F
        end

        # grids = kF_densed_kgrids(dispersion=dispersion, anglemesh=theta,
        #    bound=[0.0, 2.0])
        # cm = CompositeMesh(theta, grids)
        DIM = 2
        lattice = Matrix([1.0 0; 0 1]')
        br = BZMeshes.Cell(lattice=lattice)

        pm = PolarMesh(dispersion=dispersion, anglemesh=theta, br=br, kmax=2.0)
        @test AbstractMeshes.volume(pm) ≈ 4π

    end

    @testset "3D" begin
        dispersion(k) = dot(k, k) - 1.0

        N = 6
        bound = [-π, π]
        phi = SimpleGrid.Uniform(bound, N; isperiodic=true)

        N = 4
        bound = [-π / 2, π / 2]
        theta = SimpleGrid.Uniform(bound, N; isperiodic=true)

        am = CompositeMesh(phi, [theta for i in 1:length(phi)])

        k_F_previous = 0.0
        for ap in am
            f(k) = dispersion(BZMeshes._spherical2cart(Spherical(k, ap...)))
            k_F = find_zero(f, k_F_previous)
            @test k_F ≈ 1.0
            @test find_kFermi(dispersion, ap; kinit=k_F_previous) ≈ 1.0
            k_F_previous = k_F
        end

        DIM = 3
        lattice = Matrix([1.0 1.0 0; 1 0 1; 0 1 1]')
        br = BZMeshes.Cell(lattice=lattice)

        pm = PolarMesh(dispersion=dispersion, anglemesh=am, br=br, kmax=2.0)
        @test AbstractMeshes.volume(pm) ≈ 32π / 3

    end

    @testset "Radial rescale" begin
        @testset "RescaledGrid" begin
            rmax = 2.0
            N = 10
            rgrid = SimpleG.Uniform([0.0, rmax^2], N; gpbound=[rmax^2 / 2N, rmax^2 * (1 - 1 / 2N)])
            rrg = radial_rescale(grid=rgrid, DIM=2)
            println(rrg.grid)
            bound = [-π, π]
            theta = SimpleGrid.Uniform(bound, N; gpbound=[π * (1 / N - 1), π * (1 - 1 / N)], isperiodic=true)

            # for (i, p) in enumerate(theta)
            #     println(volume(theta, i))
            # end

            DIM = 2
            lattice = Matrix([1.0 0; 0 1]')
            br = BZMeshes.Brillouin(lattice=lattice)
            bzmesh = PolarMesh(br, CompositeMesh(theta, [rrg for i in 1:length(theta)]))
            for (i, p) in enumerate(bzmesh)
                println(p, AbstractMeshes.volume(bzmesh, i))
            end
        end

    end
end
