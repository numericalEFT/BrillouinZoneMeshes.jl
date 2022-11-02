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

@testset "PolarMesh Generator" begin
    @testset "Find k_F" begin
        # given dispersion function accept a k in cartesian
        # goal is to find k_F at direction specified by angle
        dispersion(k) = dot(k, k) - 1.0

        # 2d
        N = 10
        bound = [0, 2π]
        theta = SimpleGrid.Uniform(bound, N; isperiodic=true)

        k_F_previous = 0.0
        for θ in theta
            f(k) = dispersion(BZMeshes._polar2cart(Polar(k, θ)))
            k_F = find_zero(f, k_F_previous)
            @test k_F ≈ 1.0
            @test find_kFermi(dispersion, θ; kinit=k_F_previous) ≈ 1.0
            k_F_previous = k_F
        end

        grids = kF_densed_kgrids(dispersion=dispersion, anglemesh=theta,
            bound=[0.0, 2.0])
        cm = CompositeMesh(theta, grids)
        DIM = 2
        lattice = Matrix([1.0 0; 0 1]')
        br = BZMeshes.Brillouin(lattice=lattice)

        pm = PolarMesh(br, cm)
        @test AbstractMeshes.volume(pm) ≈ 4π
    end
end
