@testset "BZMeshes.jl" begin
    using BrillouinZoneMeshes.BZMeshes
    using BrillouinZoneMeshes.BZMeshes.Coordinates
    using BrillouinZoneMeshes.CompositeGrids
    using BrillouinZoneMeshes.BaseMesh
    using BrillouinZoneMeshes.AbstractMeshes
    using BrillouinZoneMeshes.LinearAlgebra
    using BrillouinZoneMeshes.StaticArrays
    using BrillouinZoneMeshes.Roots
    using BrillouinZoneMeshes.BZMeshes: find_kFermi, find_kFermi, radial_rescale, kF_densed_kgrids
    rng = MersenneTwister(1234)

    @testset "UniformBZMesh" begin
        DIM = 2
        N1, N2 = 4, 5
        lattice = Matrix([1.0 0; 0 1]')
        br = BZMeshes.Cell(lattice=lattice)
        mesh = UniformBZMesh(cell=br, size=(N1, N2))

        @test length(mesh) == N1 * N2
        @test size(mesh) == (N1, N2)

        vol = 0.0
        for (i, x) in enumerate(mesh)
            fracx = mesh[AbstractMeshes.FracCoords, i]
            @test fracx ≈ AbstractMeshes.cart_to_frac(mesh, x)
            @test AbstractMeshes.locate(mesh, x) == i
            vol += AbstractMeshes.volume(mesh, i)
        end
        @test AbstractMeshes.volume(mesh) ≈ vol
    end

    @testset "PolarMeshes" begin
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
                cm = ProdMesh(grids, theta)

                DIM = 2
                lattice = Matrix([1.0 0; 0 1]')
                br = BZMeshes.Cell(lattice=lattice)

                pm = PolarMesh(br, cm)
                vol = 0.0
                for (i, p) in enumerate(pm)
                    @test p == BZMeshes._polar2cart(pm[AngularCoords, i])
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
                am = ProdMesh([theta for i in 1:length(phi)], phi)
                println(typeof(size(am)))
                cm = ProdMesh([rg for i in 1:length(am)], am)
                println(typeof(size(cm)))
                println(typeof(size(cm.mesh)))

                DIM = 3
                lattice = Matrix([1.0 0 0; 0 1 0; 0 0 1]')
                br = BZMeshes.Cell(lattice=lattice)

                pm = PolarMesh(br, cm)
                println(typeof(size(pm)))

                vol = 0.0
                for (i, p) in enumerate(pm)
                    @test p == BZMeshes._spherical2cart(pm[AngularCoords, i])
                    @test AbstractMeshes.locate(pm, p) == i
                    vol += AbstractMeshes.volume(pm, i)
                end
                # a ball with r=2.0
                @test vol ≈ 32π / 3
            end
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
                # cm = ProdMesh(theta, grids)
                DIM = 2
                lattice = Matrix([1.0 0; 0 1]')
                br = BZMeshes.Cell(lattice=lattice)

                pm = PolarMesh(dispersion=dispersion, anglemesh=theta, cell=br, kmax=2.0)
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

                am = ProdMesh(phi, [theta for i in 1:length(phi)])

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

                pm = PolarMesh(dispersion=dispersion, anglemesh=am, cell=br, kmax=2.0)
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
                    br = BZMeshes.Cell(lattice=lattice)
                    bzmesh = PolarMesh(br, ProdMesh(theta, [rrg for i in 1:length(theta)]))
                    for (i, p) in enumerate(bzmesh)
                        println(p, AbstractMeshes.volume(bzmesh, i))
                    end
                end

            end
        end

    end
end

