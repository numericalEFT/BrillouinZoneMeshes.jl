@testset "PolarMeshes" begin
    using BrillouinZoneMeshes.BZMeshes
    using BrillouinZoneMeshes.BZMeshes.Coordinates
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
                [-π, π],
                [-π, 0.0, π],
                N,
                0.1,
                M
            )
            println(theta)
            grids = [CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.0], [sqrt(a * cos(θ)^2 + b * sin(θ)^2),], N, 0.1, M) for θ in theta]
            cm = CompositeMesh(theta, grids)

            DIM = 2
            lattice = Matrix([1.0 0; 0 1]')
            br = BZMeshes.Brillouin(lattice=lattice)

            pm = PolarMesh(br, cm)
            vol = 0.0
            for (i, p) in enumerate(pm)
                @test p == BZMeshes._polar2cart(pm[Angular, i])
                @test AbstractMeshes.locate(pm, p) == i
                vol += AbstractMeshes.volume(pm, i)
            end
            # grid is on a circle with r=2.0
            @test vol ≈ 4π
        end

        @testset "3D PolarMesh" begin
            N, M = 2, 2
            # theta grid dense around 0 and π
            phi = CompositeGrid.LogDensedGrid(
                :cheb,
                [-π, π],
                [-π, 0.0, π],
                N,
                0.1,
                M
            )
            theta = CompositeGrid.LogDensedGrid(
                :cheb,
                [-π / 2, π / 2],
                [0.0,],
                N,
                0.1,
                M
            )
            rg = CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.0], [1.0,], N, 0.1, M)
            am = CompositeMesh(phi, [theta for i in 1:length(phi)])
            println(typeof(size(am)))
            cm = CompositeMesh(am, [rg for i in 1:length(am)])
            println(typeof(size(cm)))
            println(typeof(size(cm.mesh)))

            DIM = 3
            lattice = Matrix([1.0 0 0; 0 1 0; 0 0 1]')
            br = BZMeshes.Brillouin(lattice=lattice)

            pm = PolarMesh(br, cm)
            println(typeof(size(pm)))

            vol = 0.0
            for (i, p) in enumerate(pm)
                @test p == BZMeshes._spherical2cart(pm[Angular, i])
                @test AbstractMeshes.locate(pm, p) == i
                vol += AbstractMeshes.volume(pm, i)
            end
            # a ball with r=2.0
            @test vol ≈ 32π / 3
        end
    end
end