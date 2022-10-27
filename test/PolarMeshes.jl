@testset "PolarMeshes" begin
    using BrillouinZoneMeshes.BZMeshes
    using BrillouinZoneMeshes.CoordinateTransformations
    using BrillouinZoneMeshes.CompositeGrids
    using BrillouinZoneMeshes.BaseMesh
    using BrillouinZoneMeshes.AbstractMeshes

    rng = MersenneTwister(1234)
    @testset "CoordinateTransformations" begin

        # special cases
        r = Polar(1, 0)
        @test BZMeshes._polar2cart(r) ≈ [1.0, 0]
        r = Polar(1, π / 2)
        @test BZMeshes._polar2cart(r) ≈ [0, 1.0]
        r = Polar(1, π)
        @test BZMeshes._polar2cart(r) ≈ [-1.0, 0]
        r = Polar(1, 3π / 2)
        @test BZMeshes._polar2cart(r) ≈ [0, -1.0]

        # random test
        Ntest = 16
        for i in 1:Ntest
            x = rand(rng, 2)
            @test BZMeshes._polar2cart(BZMeshes._cart2polar(x)) ≈ x
        end

        Ntest = 16
        for i in 1:Ntest
            x = rand(rng, 3)
            @test BZMeshes._spherical2cart(BZMeshes._cart2spherical(x)) ≈ x
        end

    end

    @testset "PolarMeshes" begin
        @testset "2D PolarMesh" begin

            a, b = 0.8, 1.2

            N, M = 3, 2
            # theta grid dense around 0 and π
            theta = CompositeGrid.LogDensedGrid(
                :cheb,
                [0.0, 2π],
                [0.0, π, 2π],
                N,
                0.1,
                M
            )
            println(theta)
            grids = [CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.0], [sqrt(a * cos(θ)^2 + b * sin(θ)^2),], N, 0.1, M) for θ in theta]
            cm = CompositeMesh(theta, grids)

        end
    end
end