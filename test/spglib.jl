
@testset "Spglib" begin
    getidx(umesh, kidx) = BrillouinZoneMeshes.spglib_grid_address_to_index(umesh, kidx)
    function test(dim, lattice, atoms, pos, ksize, kshift)
        println("Test $dim-D grid size: ", ksize, " with a shift ", kshift)
        _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]
        _kshift = _kshift[1:dim]
        br = BrillouinZoneMeshes.Model.Brillouin(lattice=lattice, atoms=atoms, positions=pos)
        brmesh = BrillouinZoneMeshes.BZMeshes.Monkhorst_Pack(br=br, size=tuple(ksize...), shift=_kshift)
        umesh = brmesh.mesh
        @test getidx(umesh, [ksize[1] - 1, 0, 0]) == getidx(umesh, [-1, 0, 0])
        @test getidx(umesh, [0, ksize[2] - 1, 0]) == getidx(umesh, [0, -1, 0])
        if dim == 3
            @test getidx(umesh, [0, 0, ksize[3] - 1]) == getidx(umesh, [0, 0, -1])
        end
        if dim == 3
            res = [getidx(umesh, [x, y, z]) for x in 0:ksize[1]-1, y in 0:ksize[2]-1, z in 0:ksize[3]-1]
            @test sort(vec(res)) == collect(1:length(umesh))
        elseif dim == 2
            res = [getidx(umesh, [x, y, 0]) for x in 0:ksize[1]-1, y in 0:ksize[2]-1]
            @test sort(vec(res)) == collect(1:length(umesh))
        end
    end
    test(3, silicon.lattice, silicon.atoms, silicon.positions, [4, 4, 4], false)
    test(3, silicon.lattice, silicon.atoms, silicon.positions, [4, 4, 4], true)

    test(3, silicon.lattice, silicon.atoms, silicon.positions, [5, 4, 3], false)
    test(3, silicon.lattice, silicon.atoms, silicon.positions, [5, 4, 3], true)

    lattice = [[1.0, 0.0] [0.0, 1.0]]

    test(2, lattice, [1,], [[0.0, 0.0],], [4, 4], false)
    test(2, lattice, [1,], [[0.0, 0.0],], [4, 4], true)

    test(2, lattice, [1,], [[0.0, 0.0],], [5, 4], false)
    test(2, lattice, [1,], [[0.0, 0.0],], [5, 4], true)
end