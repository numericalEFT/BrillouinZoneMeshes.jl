@testset "Base Mesh" begin
    rng = MersenneTwister(1234)

    @testset "UMesh" begin
        DIM = 2
        N1, N2 = 3, 5
        lattice = Matrix([1/N1/2 0; 0 1.0/N2/2]') .* 2π
        # so that bzmesh[i,j] = (2i-1,2j-1)
        cell = BZMeshes.Cell(lattice=lattice)
        mesh = BaseMesh.UMesh(br=cell, origin=ones(DIM) ./ 2, size=(N1, N2), shift=zeros(DIM))

        @inferred mesh[1]
        @inferred mesh[1, 1]
        @inferred mesh[AbstractMeshes.FracCoords, 1]
        @inferred AbstractMeshes.locate(mesh, [0.0, 0.0])
        @inferred AbstractMeshes.volume(mesh, 1)
        @inferred AbstractMeshes.volume(mesh)

        @test length(mesh) == N1 * N2
        @test size(mesh) == (N1, N2)

        vol = 0.0
        for (i, x) in enumerate(mesh)
            fracx = mesh[AbstractMeshes.FracCoords, i]
            @test mesh[i] ≈ x
            @test fracx ≈ AbstractMeshes.cart_to_frac(mesh, x)
            @test AbstractMeshes.locate(mesh, x) == i
            vol += AbstractMeshes.volume(mesh, i)
        end
        @test AbstractMeshes.volume(mesh) ≈ vol
        data = ones(Float64, size(mesh))
        @test AbstractMeshes.interp(data, mesh, [0.0, 0.0]) ≈ 1.0
        @test AbstractMeshes.integrate(data, mesh) ≈ vol
        @test AbstractMeshes.lattice_vector(mesh, 1) ≈ [2 * N1, 0]
        @test AbstractMeshes.inv_lattice_vector(mesh, 1) ≈ [1 / 2 / N1, 0]
    end

    @testset "ProdMesh" begin
        using BrillouinZoneMeshes.CompositeGrids
        using BrillouinZoneMeshes.BaseMesh
        using BrillouinZoneMeshes.AbstractMeshes

        # basics
        struct NotAPM{T,DIM} <: BaseMesh.AbstractProdMesh{T,DIM} end
        @test_throws ErrorException BaseMesh._getgrid(NotAPM{Int,3}(), 1)

        @testset "DirectProdMesh" begin
            N, M = 3, 2
            r = CompositeGrid.LogDensedGrid(
                :cheb,
                [0, 2.0],
                [1.0,],
                N,
                0.01,
                M
            )
            phi = SimpleGrid.Uniform([-π, π], 6)
            theta = SimpleGrid.Uniform([-π / 2, π / 2], 4)
            dpm = DirectProdMesh(r, theta, phi)
            println(size(dpm))

            @inferred dpm[1]

            vol = 0.0
            for (pi, p) in enumerate(dpm)
                i, j, k = AbstractMeshes._ind2inds(size(dpm), pi)
                @test dpm[pi] ≈ [r[i], theta[j], phi[k]]
                @test pi == AbstractMeshes.locate(dpm, p)
                vol += AbstractMeshes.volume(dpm, pi)
            end
            @test vol ≈ AbstractMeshes.volume(dpm)
        end
        @testset "Construct ProdMesh" begin

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

            cm = ProdMesh(grids, theta)
            println([cm.grids[i].panel[2] for i in 1:length(theta)])
            println(size(cm))

            @inferred cm[1]

            for j in 1:length(cm.mesh)
                for i in 1:length(cm.grids[j])
                    p = cm[i, j]
                    @test p[1] == cm.grids[j][i]
                    @test p[2] == cm.mesh[j]
                end
            end
            vol = 0.0
            for (pi, p) in enumerate(cm)
                # println(pi, p)
                @test pi == AbstractMeshes.locate(cm, p)
                vol += AbstractMeshes.volume(cm, pi)
            end
            @test vol ≈ AbstractMeshes.volume(cm)
        end
    end

    @testset "ChebMesh" begin
        # 2d
        N, DIM = 10, 2
        origin = [0.0 0.0]
        latvec = [π 0; 0 π]'
        cm = ChebMesh(origin, latvec, DIM, N)

        cm2 = ChebMesh(origin, latvec, cm)

        @inferred cm[1]

        vol = 0.0
        for (i, x) in enumerate(cm2)
            @test cm2[i] ≈ x
            @test AbstractMeshes.locate(cm2, x) == i
            vol += AbstractMeshes.volume(cm2, i)
        end
        @test AbstractMeshes.volume(cm2) ≈ vol
        # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
        f(x) = sin(x[1]) + cos(x[2])

        data = zeros(Float64, size(cm))

        for i in 1:length(cm)
            data[i] = f(cm[i])
        end

        ## interpolate
        testN = 3
        xlist = rand(rng, testN) * π
        ylist = rand(rng, testN) * π
        for x in xlist
            for y in ylist
                @test isapprox(f([x, y]), interp(data, cm, [x, y]), rtol=4e-2)
            end
        end

        ## integrate
        integral = integrate(data, cm)
        # println("integral=$(integral)")
        @test isapprox(integral, 2π, rtol=1e-3)

        # 3d
        N, DIM = 10, 3
        origin = [0.0 0.0 0.0]
        latvec = [1 0 0; 0 1 0; 0 0 1]'
        cm = ChebMesh(origin, latvec, DIM, N)

        g(x) = exp(-dot(x, x))

        data = zeros(Float64, size(cm))

        for i in 1:length(cm)
            data[i] = g(cm[i])
        end
        ## interpolate
        testN = 3
        xlist = rand(rng, testN)
        ylist = rand(rng, testN)
        zlist = rand(rng, testN)
        for x in xlist
            for y in ylist
                for z in zlist
                    @test isapprox(g([x, y, z]), interp(data, cm, [x, y, z]), rtol=4e-2)
                end
            end
        end
    end

    ##########################
    # LEGACY BELOW
    #############################

    # function dispersion(k)
    #     me = 0.5
    #     return dot(k, k) / 2me
    # end

    # function density(k)
    #     T = 0.01
    #     μ = 1.0

    #     ϵ = dispersion(k) - μ

    #     # return 1 / (exp((ϵ) / T) + 1.0)
    #     return (π * T)^2 / ((π * T)^2 + ϵ^2)
    #     # return (exp((ϵ) / T) / T)/(exp((ϵ) / T) + 1.0)^2
    # end

    # @testset "Centered Uniform Mesh" begin
    #     δ = 1e-3

    #     # index convert

    #     # 2D
    #     N, DIM = 4, 2
    #     origin = [0.0 0.0]
    #     latvec = [2 0; 1 sqrt(3)]'
    #     # latvec = [1 0; 0 1]'

    #     umesh = UniformMesh{DIM,N}(origin, latvec)
    #     for i in 1:N^(DIM)
    #         # println("$i -> $(BaseMesh._ind2inds(i, N, DIM))")
    #         @test i == BaseMesh._inds2ind(umesh, BaseMesh._ind2inds(umesh, i))
    #     end

    #     for i in 1:length(umesh)
    #         @test umesh[i] == umesh[BaseMesh._ind2inds(umesh, i)...]
    #     end

    #     for i in 1:length(umesh)
    #         for j in 1:DIM
    #             shift = zeros(Float64, DIM)
    #             indshift = zeros(Int, DIM)

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             # println(inds)

    #             shift = δ .* latvec[:, j]
    #             indshift[j] = 0
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             shift = -δ .* latvec[:, j]
    #             indshift[j] = -1
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             if inds[j] == 1
    #                 indshift[j] = 0
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))
    #         end
    #     end

    #     # 3D
    #     N, DIM = 3, 3
    #     origin = [0.0 0.0 0.0]
    #     latvec = [1.0 0 0; 0 1.0 0; 0 0 1.0]'

    #     umesh = UniformMesh{DIM,N}(origin, latvec)

    #     for i in 1:length(umesh)
    #         for j in 1:DIM
    #             shift = zeros(Float64, DIM)
    #             indshift = zeros(Int, DIM)

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             # println(inds)

    #             shift = δ .* latvec[:, j]
    #             indshift[j] = 0
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             shift = -δ .* latvec[:, j]
    #             indshift[j] = -1
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             if inds[j] == 1
    #                 indshift[j] = 0
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))
    #         end
    #     end

    # end

    # @testset "Edged Uniform Mesh" begin
    #     δ = 1e-3

    #     # 2D
    #     N, DIM = 4, 2
    #     origin = [0.0 0.0]
    #     latvec = [2 0; 1 sqrt(3)]'
    #     # latvec = [1 0; 0 1]'

    #     umesh = UniformMesh{DIM,N,BaseMesh.EdgedMesh}(origin, latvec)

    #     for i in 1:length(umesh)
    #         for j in 1:DIM
    #             shift = zeros(Float64, DIM)
    #             indshift = zeros(Int, DIM)

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             # println(inds)

    #             shift = δ .* latvec[:, j]
    #             indshift[j] = 0
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             # @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)

    #             inds = BaseMesh._ind2inds(umesh, i)
    #             inds = Vector(inds)
    #             shift = -δ .* latvec[:, j]
    #             indshift[j] = -1
    #             for k in 1:DIM
    #                 if inds[k] == N
    #                     inds[k] = N - 1
    #                     if k == j
    #                         indshift[j] = 0
    #                     end
    #                 end
    #             end
    #             if inds[j] == 1
    #                 indshift[j] = 0
    #             end
    #             ind = BaseMesh._inds2ind(umesh, inds + indshift)
    #             # @test ind == floor(umesh, umesh[i] + shift)
    #             # @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)
    #         end
    #     end

    # end

    # @testset "Interp and Integral for Uniform" begin
    #     # 2d
    #     N, DIM = 100, 2
    #     origin = [0.0 0.0]
    #     latvec = [π 0; 0 π]'
    #     umesh = UniformMesh{DIM,N}(origin, latvec)

    #     # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
    #     f(x) = sin(x[1]) + cos(x[2])

    #     data = zeros(Float64, size(umesh))

    #     for i in 1:length(umesh)
    #         data[i] = f(umesh[i])
    #     end

    #     ## interpolate
    #     testN = 3
    #     xlist = rand(rng, testN) * π
    #     ylist = rand(rng, testN) * π
    #     for x in xlist
    #         for y in ylist
    #             @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=9e-2)
    #         end
    #     end

    #     ## integrate
    #     integral = BaseMesh.integrate(data, umesh)
    #     # println("integral=$(integral)")
    #     @test isapprox(integral, 2π, rtol=1e-3)

    #     # 2d edgedmesh
    #     N, DIM = 100, 2
    #     origin = [0.0 0.0]
    #     latvec = [π 0; 0 π]'
    #     umesh = UniformMesh{DIM,N,BaseMesh.EdgedMesh}(origin, latvec)

    #     data = zeros(Float64, size(umesh))

    #     for i in 1:length(umesh)
    #         data[i] = f(umesh[i])
    #     end

    #     ## interpolate
    #     testN = 3
    #     xlist = rand(rng, testN) * π
    #     ylist = rand(rng, testN) * π
    #     for x in xlist
    #         for y in ylist
    #             @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=9e-2)
    #         end
    #     end

    #     ## integrate
    #     integral = BaseMesh.integrate(data, umesh)
    #     # println("integral=$(integral)")
    #     @test isapprox(integral, 2π, rtol=1e-3)

    #     # 3d
    #     N, DIM = 100, 3
    #     origin = [0.0 0.0 0.0]
    #     latvec = [1 0 0; 0 1 0; 0 0 1]'
    #     umesh = UniformMesh{DIM,N}(origin, latvec)

    #     g(x) = x[1] * x[2] * x[3]

    #     data = zeros(Float64, size(umesh))

    #     for i in 1:length(umesh)
    #         data[i] = f(umesh[i])
    #     end
    #     ## interpolate
    #     testN = 3
    #     xlist = rand(rng, testN)
    #     ylist = rand(rng, testN)
    #     zlist = rand(rng, testN)
    #     for x in xlist
    #         for y in ylist
    #             for z in zlist
    #                 @test isapprox(f([x, y, z]), BaseMesh.interp(data, umesh, [x, y, z]), rtol=4e-2)
    #             end
    #         end
    #     end

    # end

    # @testset "Interp and Integral for BaryCheb" begin
    #     # 2d
    #     N, DIM = 20, 2
    #     origin = [0.0 0.0]
    #     latvec = [π 0; 0 π]'
    #     umesh = BaryChebMesh(origin, latvec, DIM, N)

    #     # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
    #     f(x) = sin(x[1]) + cos(x[2])

    #     data = zeros(Float64, size(umesh))

    #     for i in 1:length(umesh)
    #         data[i] = f(umesh[i])
    #     end

    #     ## interpolate
    #     testN = 3
    #     xlist = rand(rng, testN) * π
    #     ylist = rand(rng, testN) * π
    #     for x in xlist
    #         for y in ylist
    #             @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=4e-2)
    #         end
    #     end

    #     ## integrate
    #     integral = BaseMesh.integrate(data, umesh)
    #     # println("integral=$(integral)")
    #     @test isapprox(integral, 2π, rtol=1e-3)

    #     # 3d
    #     N, DIM = 100, 3
    #     origin = [0.0 0.0 0.0]
    #     latvec = [1 0 0; 0 1 0; 0 0 1]'
    #     umesh = UniformMesh{DIM,N}(origin, latvec)

    #     g(x) = x[1] * x[2] * x[3]

    #     data = zeros(Float64, size(umesh))

    #     for i in 1:length(umesh)
    #         data[i] = f(umesh[i])
    #     end
    #     ## interpolate
    #     testN = 3
    #     xlist = rand(rng, testN)
    #     ylist = rand(rng, testN)
    #     zlist = rand(rng, testN)
    #     for x in xlist
    #         for y in ylist
    #             for z in zlist
    #                 @test isapprox(f([x, y, z]), BaseMesh.interp(data, umesh, [x, y, z]), rtol=4e-2)
    #             end
    #         end
    #     end

    # end

end
