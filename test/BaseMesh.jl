@testset "Base Mesh" begin
    rng = MersenneTwister(1234)

    @testset "Brillouin" begin

        # square lattice
        DIM = 2
        lattice = Matrix([1.0 0; 0 1]')
        br = BaseMesh.Brillouin(lattice=lattice)
        @test br.inv_lattice .* 2π ≈ br.recip_lattice
        @test br.unit_cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM

        # triagular lattice
        DIM = 2
        lattice = Matrix([2.0 0; 1 sqrt(3)]')
        br = BaseMesh.Brillouin(lattice=lattice)
        @test br.inv_lattice .* 2π ≈ br.recip_lattice
        @test br.unit_cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM

        # 3d testing lattice
        DIM = 3
        lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
        br = BaseMesh.Brillouin(lattice=lattice)
        @test br.inv_lattice .* 2π ≈ br.recip_lattice
        @test br.unit_cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM
    end

    @testset "UniformBZMesh" begin
        @testset "Indexing" begin
            size = (3, 4, 5)
            for i in 1:prod(size)
                @test i == BZMeshes._inds2ind(size, BZMeshes._ind2inds(size, i))
            end
        end

        @testset "Array Interface" begin
            N1, N2 = 3, 5
            lattice = Matrix([1/N1/2 0; 0 1.0/N2/2]') .* 2π
            # so that bzmesh[i,j] = (2i-1,2j-1)
            br = BZMeshes.Brillouin(lattice=lattice)
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=(N1, N2), origin=0)
            for (pi, p) in enumerate(bzmesh)
                @test bzmesh[pi] ≈ p # linear index
                inds = BZMeshes._ind2inds(size(bzmesh), pi)
                @test p ≈ inds .* 2.0 .- 1.0
                @test bzmesh[inds...] ≈ p # cartesian index
            end
        end

        @testset "locate and volume" begin
            msize = (3, 5, 7)
            lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
            # size = (3, 5)
            # lattice = Matrix([2.3 0; 0 7.0]')
            br = BZMeshes.Brillouin(lattice=lattice)
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=msize)
            vol = 0.0
            for (pi, p) in enumerate(bzmesh)
                @test bzmesh[pi] ≈ p # linear index
                inds = BZMeshes._ind2inds(size(bzmesh), pi)
                @test bzmesh[inds...] ≈ p # cartesian index

                @test AbstractMeshes.locate(bzmesh, p) == pi
                vol += AbstractMeshes.volume(bzmesh, pi)
            end
            @test vol ≈ AbstractMeshes.volume(bzmesh)
        end

        @testset "origin and shift convention" begin
            DIM = 2
            lattice = Matrix([1.0 0; 0 1]')
            br = BZMeshes.Brillouin(lattice=lattice)

            # even numbers
            size = (4, 4)
            # Gamma-centered, no shift
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=size, origin=0, shift=0)
            @test bzmesh[1, 1] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])
            # M-P, no shift
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=size, origin=-1 / 2, shift=0)
            @test bzmesh[3, 3] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])

            # odd numbers
            size = (5, 5)
            # Gamma-centered, no shift
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=size, origin=0, shift=0)
            @test bzmesh[1, 1] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])
            # M-P, 1/2 shift
            bzmesh = BZMeshes.UniformBZMesh(br=br, size=size, origin=-1 / 2, shift=1 // 2)
            @test bzmesh[3, 3] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])

        end
    end


    function dispersion(k)
        me = 0.5
        return dot(k, k) / 2me
    end

    function density(k)
        T = 0.01
        μ = 1.0

        ϵ = dispersion(k) - μ

        # return 1 / (exp((ϵ) / T) + 1.0)
        return (π * T)^2 / ((π * T)^2 + ϵ^2)
        # return (exp((ϵ) / T) / T)/(exp((ϵ) / T) + 1.0)^2
    end

    @testset "Centered Uniform Mesh" begin
        δ = 1e-3

        # index convert

        # 2D
        N, DIM = 4, 2
        origin = [0.0 0.0]
        latvec = [2 0; 1 sqrt(3)]'
        # latvec = [1 0; 0 1]'

        umesh = UniformMesh{DIM,N}(origin, latvec)
        for i in 1:N^(DIM)
            # println("$i -> $(BaseMesh._ind2inds(i, N, DIM))")
            @test i == BaseMesh._inds2ind(umesh, BaseMesh._ind2inds(umesh, i))
        end

        for i in 1:length(umesh)
            @test umesh[i] == umesh[BaseMesh._ind2inds(umesh, i)...]
        end

        for i in 1:length(umesh)
            for j in 1:DIM
                shift = zeros(Float64, DIM)
                indshift = zeros(Int, DIM)

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                # println(inds)

                shift = δ .* latvec[:, j]
                indshift[j] = 0
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                shift = -δ .* latvec[:, j]
                indshift[j] = -1
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                if inds[j] == 1
                    indshift[j] = 0
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))
            end
        end

        # 3D
        N, DIM = 3, 3
        origin = [0.0 0.0 0.0]
        latvec = [1.0 0 0; 0 1.0 0; 0 0 1.0]'

        umesh = UniformMesh{DIM,N}(origin, latvec)

        for i in 1:length(umesh)
            for j in 1:DIM
                shift = zeros(Float64, DIM)
                indshift = zeros(Int, DIM)

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                # println(inds)

                shift = δ .* latvec[:, j]
                indshift[j] = 0
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                shift = -δ .* latvec[:, j]
                indshift[j] = -1
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                if inds[j] == 1
                    indshift[j] = 0
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                @test BaseMesh._ind2inds(umesh, ind) == BaseMesh._ind2inds(umesh, floor(umesh, umesh[i] + shift))
            end
        end

    end

    @testset "Edged Uniform Mesh" begin
        δ = 1e-3

        # 2D
        N, DIM = 4, 2
        origin = [0.0 0.0]
        latvec = [2 0; 1 sqrt(3)]'
        # latvec = [1 0; 0 1]'

        umesh = UniformMesh{DIM,N,BaseMesh.EdgedMesh}(origin, latvec)

        for i in 1:length(umesh)
            for j in 1:DIM
                shift = zeros(Float64, DIM)
                indshift = zeros(Int, DIM)

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                # println(inds)

                shift = δ .* latvec[:, j]
                indshift[j] = 0
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                # @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)

                inds = BaseMesh._ind2inds(umesh, i)
                inds = Vector(inds)
                shift = -δ .* latvec[:, j]
                indshift[j] = -1
                for k in 1:DIM
                    if inds[k] == N
                        inds[k] = N - 1
                        if k == j
                            indshift[j] = 0
                        end
                    end
                end
                if inds[j] == 1
                    indshift[j] = 0
                end
                ind = BaseMesh._inds2ind(umesh, inds + indshift)
                # @test ind == floor(umesh, umesh[i] + shift)
                # @test BaseMesh._ind2inds(ind, N, DIM) == BaseMesh._ind2inds(floor(umesh, umesh[i] + shift), N, DIM)
            end
        end

    end

    @testset "Interp and Integral for Uniform" begin
        # 2d
        N, DIM = 100, 2
        origin = [0.0 0.0]
        latvec = [π 0; 0 π]'
        umesh = UniformMesh{DIM,N}(origin, latvec)

        # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
        f(x) = sin(x[1]) + cos(x[2])

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end

        ## interpolate
        testN = 3
        xlist = rand(rng, testN) * π
        ylist = rand(rng, testN) * π
        for x in xlist
            for y in ylist
                @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=9e-2)
            end
        end

        ## integrate
        integral = BaseMesh.integrate(data, umesh)
        # println("integral=$(integral)")
        @test isapprox(integral, 2π, rtol=1e-3)

        # 2d edgedmesh
        N, DIM = 100, 2
        origin = [0.0 0.0]
        latvec = [π 0; 0 π]'
        umesh = UniformMesh{DIM,N,BaseMesh.EdgedMesh}(origin, latvec)

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end

        ## interpolate
        testN = 3
        xlist = rand(rng, testN) * π
        ylist = rand(rng, testN) * π
        for x in xlist
            for y in ylist
                @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=9e-2)
            end
        end

        ## integrate
        integral = BaseMesh.integrate(data, umesh)
        # println("integral=$(integral)")
        @test isapprox(integral, 2π, rtol=1e-3)

        # 3d
        N, DIM = 100, 3
        origin = [0.0 0.0 0.0]
        latvec = [1 0 0; 0 1 0; 0 0 1]'
        umesh = UniformMesh{DIM,N}(origin, latvec)

        g(x) = x[1] * x[2] * x[3]

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end
        ## interpolate
        testN = 3
        xlist = rand(rng, testN)
        ylist = rand(rng, testN)
        zlist = rand(rng, testN)
        for x in xlist
            for y in ylist
                for z in zlist
                    @test isapprox(f([x, y, z]), BaseMesh.interp(data, umesh, [x, y, z]), rtol=4e-2)
                end
            end
        end

    end

    @testset "Interp and Integral for BaryCheb" begin
        # 2d
        N, DIM = 20, 2
        origin = [0.0 0.0]
        latvec = [π 0; 0 π]'
        umesh = BaryChebMesh(origin, latvec, DIM, N)

        # f(x) = x[1] + 2 * x[2] + x[1] * x[2]
        f(x) = sin(x[1]) + cos(x[2])

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end

        ## interpolate
        testN = 3
        xlist = rand(rng, testN) * π
        ylist = rand(rng, testN) * π
        for x in xlist
            for y in ylist
                @test isapprox(f([x, y]), BaseMesh.interp(data, umesh, [x, y]), rtol=4e-2)
            end
        end

        ## integrate
        integral = BaseMesh.integrate(data, umesh)
        # println("integral=$(integral)")
        @test isapprox(integral, 2π, rtol=1e-3)

        # 3d
        N, DIM = 100, 3
        origin = [0.0 0.0 0.0]
        latvec = [1 0 0; 0 1 0; 0 0 1]'
        umesh = UniformMesh{DIM,N}(origin, latvec)

        g(x) = x[1] * x[2] * x[3]

        data = zeros(Float64, size(umesh))

        for i in 1:length(umesh)
            data[i] = f(umesh[i])
        end
        ## interpolate
        testN = 3
        xlist = rand(rng, testN)
        ylist = rand(rng, testN)
        zlist = rand(rng, testN)
        for x in xlist
            for y in ylist
                for z in zlist
                    @test isapprox(f([x, y, z]), BaseMesh.interp(data, umesh, [x, y, z]), rtol=4e-2)
                end
            end
        end

    end

end
