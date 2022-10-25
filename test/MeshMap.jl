@testset "MeshMaps" begin
    @testset "MeshMap" begin

        # test MeshMap constructor
        map = [1, 2, 2, 1, 2, 6, 6, 2, 2, 6, 6, 2, 1, 2, 2, 1]
        mm = MeshMap(map)
        # println(mm.irreducible_indices)
        @test mm.irreducible_indices == [1, 2, 6]
        @test mm.inv_map[1] == [1, 4, 13, 16]
        @test mm.inv_map[2] == [2, 3, 5, 8, 9, 12, 14, 15]
        @test mm.inv_map[6] == [6, 7, 10, 11]
    end
end