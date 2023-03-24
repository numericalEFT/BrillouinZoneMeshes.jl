@testset "BZMeshes.jl" begin
    using BrillouinZoneMeshes.BZMeshes
    # using BrillouinZoneMeshes.BZMeshes.Coordinates
    using BrillouinZoneMeshes.CompositeGrids
    using BrillouinZoneMeshes.BaseMesh
    using BrillouinZoneMeshes.AbstractMeshes
    using BrillouinZoneMeshes.LinearAlgebra
    using BrillouinZoneMeshes.StaticArrays
    using BrillouinZoneMeshes.Roots
    # using BrillouinZoneMeshes.BZMeshes: find_kFermi, find_kFermi, radial_rescale, kF_densed_kgrids
    rng = MersenneTwister(1234)

    @testset "UniformBZMesh" begin
        @testset "Indexing" begin
            size = (3, 4, 5)
            for i in 1:prod(size)
                @test i == AbstractMeshes._inds2ind(size, AbstractMeshes._ind2inds(size, i))
            end
        end

        @testset "Array Interface" begin
            N1, N2 = 3, 5
            lattice = Matrix([1/N1/2 0; 0 1.0/N2/2]') .* 2π
            # so that bzmesh[i,j] = (2i-1,2j-1)
            br = BZMeshes.Cell(lattice=lattice)
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=(N1, N2), origin=0)

            println(bzmesh)
            display(bzmesh)

            for (pi, p) in enumerate(bzmesh)
                @test bzmesh[pi] ≈ p # linear index
                inds = AbstractMeshes._ind2inds(size(bzmesh), pi)
                @test p ≈ inds .* 2.0 .- 1.0
                @test bzmesh[inds...] ≈ p # cartesian index
            end
        end

        @testset "lattice vector and inverse lattice vector" begin
            lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
            msize = (3, 5, 7)
            br = BZMeshes.Cell(lattice=lattice)
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=msize)

            @test lattice_vector(bzmesh, 1) ≈ br.recip_lattice[:, 1]
            @test lattice_vector(bzmesh, 2) ≈ br.recip_lattice[:, 2]
            @test lattice_vector(bzmesh, 3) ≈ br.recip_lattice[:, 3]

            @test inv_lattice_vector(bzmesh, 1) ≈ br.inv_recip_lattice[:, 1]
            @test inv_lattice_vector(bzmesh, 2) ≈ br.inv_recip_lattice[:, 2]
            @test inv_lattice_vector(bzmesh, 3) ≈ br.inv_recip_lattice[:, 3]

            @test lattice_vector(bzmesh) * inv_lattice_vector(bzmesh) ≈ Matrix(I, 3, 3)
        end

        @testset "locate and volume" begin
            msize = (3, 5, 7)
            lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
            # size = (3, 5)
            # lattice = Matrix([2.3 0; 0 7.0]')
            br = BZMeshes.Cell(lattice=lattice)
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=msize)
            vol = 0.0
            for (i, p) in enumerate(bzmesh)
                @test bzmesh[i] ≈ p # linear index
                inds = AbstractMeshes._ind2inds(size(bzmesh), i)
                @test bzmesh[inds...] ≈ p # cartesian index

                @test AbstractMeshes.locate(bzmesh, p) == i
                vol += AbstractMeshes.volume(bzmesh, i)
            end
            @test vol ≈ AbstractMeshes.volume(bzmesh)
        end

        @testset "origin and shift convention" begin
            DIM = 2
            lattice = Matrix([1.0 0; 0 1]')
            br = BZMeshes.Cell(lattice=lattice)

            # even numbers
            size = (4, 4)
            # Gamma-centered, no shift
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=size, origin=0, shift=[false, false])
            @test bzmesh[1, 1] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])
            # M-P, no shift
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=size, origin=-1 / 2, shift=[false, false])
            @test bzmesh[3, 3] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])

            # odd numbers
            size = (5, 5)
            # Gamma-centered, no shift
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=size, origin=0, shift=[false, false])
            @test bzmesh[1, 1] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])
            # M-P, 1/2 shift
            bzmesh = BZMeshes.UniformBZMesh(cell=br, size=size, origin=-1 / 2, shift=[true, true])
            @test bzmesh[3, 3] ≈ BZMeshes.SVector{DIM,eltype(lattice)}([0.0, 0.0])

        end
    end

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

end

