using BrillouinZoneMeshes
using BrillouinZoneMeshes.BZMeshes
using BrillouinZoneMeshes.BZMeshes.Coordinates
using BrillouinZoneMeshes.CompositeGrids
using BrillouinZoneMeshes.BaseMesh
using BrillouinZoneMeshes.AbstractMeshes
using BrillouinZoneMeshes.LinearAlgebra
using BrillouinZoneMeshes.StaticArrays

using Test

function AbstractMeshes.interval(mesh::ProdMesh, I::Int)
    inds = AbstractMeshes._ind2inds(size(mesh), I)
    J = Base._sub2ind(size(mesh)[2:end], inds[2:end]...)
    x1, x2 = AbstractMeshes.interval(mesh.grids[J], inds[1])
    if length(inds) > 2
        return [(x1, x2), AbstractMeshes.interval(mesh.mesh, J)...]
    else
        return [(x1, x2), AbstractMeshes.interval(mesh.mesh, J)]
    end
end

struct CompositeMesh{T,PM,SM} <: AbstractMesh{T,1}
    panelmesh::PM
    submeshes::Vector{SM}
end

function CompositeMesh(panelmesh::PM, N) where {PM}
    T = eltype(panelmesh)
    submeshes = []
    for (i, p) in enumerate(panelmesh)
        intervals = AbstractMeshes.interval(panelmesh, i)
        DIM = length(intervals)
        origin = [intervals[j][1] for j in 1:DIM]
        lattice = diagm(DIM, DIM, [intervals[j][2] - intervals[j][1] for j in 1:DIM])
        cm = ChebMesh(origin, lattice, DIM, N)
        push!(submeshes, cm)
    end

    SM = typeof(submeshes[1])
    return CompositeMesh{T,PM,SM}(panelmesh, submeshes)
end


@testset "CompositeMesh" begin
    @testset "Constructor" begin
        a, b = 0.8, 1.2

        N, M = 3, 2
        # theta grid dense around 0 and π
        theta = CompositeGrid.LogDensedGrid(
            :cheb,
            [-π, π],
            [-π, 0, π],
            N,
            0.1,
            M
        )
        println(theta)
        grids = [CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.0], [sqrt(a * cos(θ)^2 + b * sin(θ)^2),], N, 0.1, M) for θ in theta]

        pm = ProdMesh(theta, grids)
        cm = CompositeMesh(pm, 3)
    end

    # @testset "2D" begin
    #     # given dispersion function accept a k in cartesian
    #     # goal is to find k_F at direction specified by angle
    #     dispersion(k) = dot(k, k) - 1.0

    #     # 2d
    #     N = 10
    #     bound = [-π, π]
    #     theta = SimpleGrid.Uniform(bound, N; isperiodic=true)

    #     k_F_previous = 0.0
    #     for θ in theta
    #         f(k) = dispersion(BZMeshes._polar2cart(Polar(k, θ)))
    #         k_F = find_zero(f, k_F_previous)
    #         @test k_F ≈ 1.0
    #         @test find_kFermi(dispersion, θ; kinit=k_F_previous) ≈ 1.0
    #         k_F_previous = k_F
    #     end

    #     # grids = kF_densed_kgrids(dispersion=dispersion, anglemesh=theta,
    #     #    bound=[0.0, 2.0])
    #     # cm = CompositeMesh(theta, grids)
    #     DIM = 2
    #     lattice = Matrix([1.0 0; 0 1]')
    #     br = BZMeshes.Cell(lattice=lattice)

    #     pm = PolarMesh(dispersion=dispersion, anglemesh=theta, cell=br, kmax=2.0)
    #     @test AbstractMeshes.volume(pm) ≈ 4π

    # end

    # @testset "3D" begin
    #     dispersion(k) = dot(k, k) - 1.0

    #     N = 6
    #     bound = [-π, π]
    #     phi = SimpleGrid.Uniform(bound, N; isperiodic=true)

    #     N = 4
    #     bound = [-π / 2, π / 2]
    #     theta = SimpleGrid.Uniform(bound, N; isperiodic=true)

    #     am = CompositeMesh(phi, [theta for i in 1:length(phi)])

    #     k_F_previous = 0.0
    #     for ap in am
    #         f(k) = dispersion(BZMeshes._spherical2cart(Spherical(k, ap...)))
    #         k_F = find_zero(f, k_F_previous)
    #         @test k_F ≈ 1.0
    #         @test find_kFermi(dispersion, ap; kinit=k_F_previous) ≈ 1.0
    #         k_F_previous = k_F
    #     end

    #     DIM = 3
    #     lattice = Matrix([1.0 1.0 0; 1 0 1; 0 1 1]')
    #     br = BZMeshes.Cell(lattice=lattice)

    #     pm = PolarMesh(dispersion=dispersion, anglemesh=am, cell=br, kmax=2.0)
    #     @test AbstractMeshes.volume(pm) ≈ 32π / 3

    # end

    # @testset "Radial rescale" begin
    #     @testset "RescaledGrid" begin
    #         rmax = 2.0
    #         N = 10
    #         rgrid = SimpleG.Uniform([0.0, rmax^2], N; gpbound=[rmax^2 / 2N, rmax^2 * (1 - 1 / 2N)])
    #         rrg = radial_rescale(grid=rgrid, DIM=2)
    #         println(rrg.grid)
    #         bound = [-π, π]
    #         theta = SimpleGrid.Uniform(bound, N; gpbound=[π * (1 / N - 1), π * (1 - 1 / N)], isperiodic=true)

    #         # for (i, p) in enumerate(theta)
    #         #     println(volume(theta, i))
    #         # end

    #         DIM = 2
    #         lattice = Matrix([1.0 0; 0 1]')
    #         br = BZMeshes.Cell(lattice=lattice)
    #         bzmesh = PolarMesh(br, CompositeMesh(theta, [rrg for i in 1:length(theta)]))
    #         for (i, p) in enumerate(bzmesh)
    #             println(p, AbstractMeshes.volume(bzmesh, i))
    #         end
    #     end

    # end
end