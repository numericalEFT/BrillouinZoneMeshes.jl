using BrillouinZoneMeshes
using BrillouinZoneMeshes.BZMeshes
using BrillouinZoneMeshes.BZMeshes.Coordinates
using BrillouinZoneMeshes.CompositeGrids
using BrillouinZoneMeshes.BaseMesh
using BrillouinZoneMeshes.AbstractMeshes
using BrillouinZoneMeshes.LinearAlgebra

using Roots
using Test

function find_kFermi(dispersion, angle; kinit=0.0)
    if length(angle) == 1
        return find_zero(k -> dispersion(BZMeshes._polar2cart(Polar(k, angle...))), kinit)
    elseif length(angle) == 2
        return find_zero(k -> dispersion(BZMeshes._spherical2cart(Spherical(k, angle...))), kinit)
    else
        error("dimension $(length(angle)+1) not implemented!")
    end
end

function kF_densed_kgrids(; dispersion,
    anglemesh,
    bound,
    basegridtype=:cheb,
    Nloggrid=3,
    minterval=0.01,
    Nbasegrid=2)
    # assume dispersion==0 has one root for each angle
    k1 = find_kFermi(dispersion, anglemesh[1])
    g1 = CompositeGrid.LogDensedGrid(basegridtype, bound, [k1,], Nloggrid, minterval, Nbasegrid)
    grids = [g1,]
    for i in 2:length(anglemesh)
        kF = find_kFermi(dispersion, anglemesh[i])
        g = CompositeGrid.LogDensedGrid(basegridtype, bound, [kF,], Nloggrid, minterval, Nbasegrid)
        push!(grids, g)
    end
    return grids
end


function BZMeshes.PolarMesh(; dispersion, anglemesh, br, kmax,
    kwargs...)
    bound = [0.0, kmax]
    grids = kF_densed_kgrids(dispersion=dispersion, anglemesh=anglemesh, bound=bound, kwargs...)
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
        br = BZMeshes.Brillouin(lattice=lattice)

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
        br = BZMeshes.Brillouin(lattice=lattice)

        pm = PolarMesh(dispersion=dispersion, anglemesh=am, br=br, kmax=2.0)
        @test AbstractMeshes.volume(pm) ≈ 32π / 3

    end
end
