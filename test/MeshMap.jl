@testset "MeshMaps" begin

    locate, volume = MeshMaps.locate, MeshMaps.volume
    @testset "MeshMap" begin

        # test MeshMap constructor
        map = [1, 2, 2, 1, 2, 6, 6, 2, 2, 6, 6, 2, 1, 2, 2, 1]
        mm = MeshMap(map)
        # println(mm.irreducible_indices)
        @test mm.irreducible_indices == [1, 2, 6]
        @test mm.inv_map[1] == [1, 4, 13, 16]
        @test mm.inv_map[2] == [2, 3, 5, 8, 9, 12, 14, 15]
        @test mm.inv_map[6] == [6, 7, 10, 11]

        @test length(mm) == 3
        @test size(mm) == (3,)
        for i in 1:length(mm.map)
            @test mm.irreducible_indices[MeshMaps.locate(mm, i)] == mm.map[i]
        end
    end

    @testset "ReducedBZMesh" begin
        DIM = 2
        lattice = Matrix([1.0 0; 0 1]')
        br = BZMeshes.Cell(lattice=lattice)
        umesh = BZMeshes.UniformBZMesh(br=br, size=(4, 4))

        # hand-made map
        map = [1, 2, 2, 1, 2, 6, 6, 2, 2, 6, 6, 2, 1, 2, 2, 1]
        mm = MeshMap(map)

        rmesh = ReducedBZMesh(umesh, mm)
        vol = 0.0
        for (i, p) in enumerate(rmesh)
            vol += AbstractMeshes.volume(rmesh, i)
        end
        @test vol ≈ AbstractMeshes.volume(rmesh)

        for (i, p) in enumerate(umesh)
            @test AbstractMeshes.locate(rmesh, p) == rmesh.meshmap.inv_indices[rmesh.meshmap[AbstractMeshes.locate(umesh, p)]]
        end
    end

end