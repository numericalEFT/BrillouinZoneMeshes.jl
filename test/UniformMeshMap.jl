@testset "Uniform MeshMap test" begin
    function test(dim, lattice, atoms, pos, ksize, kshift; kcoords, kweights=nothing)
        println("Test $dim-D grid size: ", ksize, " with a shift ", kshift)
        _kshift = [kshift, kshift, kshift]
        _kshift = _kshift[1:dim]
        br = BrillouinZoneMeshes.Cells.Cell(lattice=lattice, atoms=atoms, positions=pos)
        brmesh = BrillouinZoneMeshes.BZMeshes.Monkhorst_Pack(br=br, size=tuple(ksize...), shift=_kshift)

        meshmap = MeshMap(brmesh)
        # println(meshmap.irreducible_indices)
        # println(kcoords)

        # for kpoint in meshmap.irreducible_indices
        #     kvec = brmesh.mesh.inv_lattice * brmesh[kpoint]
        #     println(kvec)
        #     # _kvec = BrillouinZoneMeshes.Cells.recip_vector_red_to_cart(br, meshmap.kcoords_global[kpoint])
        #     # @test kvec ≈ _kvec
        # end
        klist = Vector{Int}()
        for coord in kcoords
            kvec = lattice_vector(brmesh) * coord
            kidx = BrillouinZoneMeshes.locate(brmesh, kvec)
            # weight = BrillouinZoneMeshes.volume(brmesh.mesh, kidx)
            # println(brmesh[kidx], " -> ", brmesh.mesh.inv_lattice * brmesh[kidx])
            push!(klist, meshmap.map[kidx])
        end
        @test sort(klist) == sort(meshmap.irreducible_indices)

        ##### test home-made irreducible_kcoord ##############
        kidxmap = PointSymmetry.irreducible_kcoord(br.symmetry, [inv_lattice_vector(brmesh) * brmesh[i] for i in 1:length(brmesh)])
        ir_klist = sort(unique(kidxmap))

        # for kpoint in ir_klist
        #     kvec = inv_lattice_vector(brmesh) * brmesh[kpoint]
        #     println(kvec)
        #     # _kvec = BrillouinZoneMeshes.Cells.recip_vector_red_to_cart(br, meshmap.kcoords_global[kpoint])
        #     # @test kvec ≈ _kvec
        # end

        @test length(ir_klist) == length(meshmap.irreducible_indices)
        for (idx, ir_idx) in enumerate(kidxmap)
            # println(idx, ": ", fractional_coordinates(brmesh, idx), " -> ", ir_idx, " : ", fractional_coordinates(brmesh, ir_idx))
            @test meshmap.map[idx] == meshmap.map[ir_idx]
        end

        return meshmap, brmesh
    end
    size, shift = [4, 4, 4], false
    kcoords = [[0.0, 0.0, 0.0],
        [0.25, 0.0, 0.0],
        [-0.5, 0.0, 0.0],
        [0.25, 0.25, 0.0],
        [-0.5, 0.25, 0.0],
        [-0.25, 0.25, 0.0],
        [-0.5, -0.5, 0.0],
        [-0.25, -0.5, 0.25]]
    kweights = [0.015625, 0.125, 0.0625, 0.09375, 0.375, 0.1875, 0.046875, 0.09375]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)
    exit(0)

    size, shift = [4, 4, 4], true
    kcoords = [
        [0.125, 0.125, 0.125],
        [0.375, 0.125, 0.125],
        [-0.375, 0.125, 0.125],
        [-0.125, 0.125, 0.125],
        [0.375, 0.375, 0.125],
        [-0.375, 0.375, 0.125],
        [-0.125, 0.375, 0.125],
        [-0.375, -0.375, 0.125],
        [0.375, 0.375, 0.375],
        [-0.375, 0.375, 0.375],
    ]
    kweights = [0.03125, 0.09375, 0.09375, 0.09375, 0.09375, 0.1875, 0.1875, 0.09375, 0.03125, 0.09375]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)

    size, shift = [3, 3, 3], false
    kcoords = [
        [0.0, 0.0, 0.0],
        [0.3333333333333333, 0.0, 0.0],
        [0.3333333333333333, 0.3333333333333333, 0.0],
        [-0.3333333333333333, 0.3333333333333333, 0.0]
    ]
    kweigths = [0.037037037037037035,
        0.2962962962962963,
        0.2222222222222222,
        0.4444444444444444]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)

    size, shift = [3, 3, 3], true
    kcoords = [
        [0.16666666666666666, 0.16666666666666666, 0.16666666666666666],
        [-0.5, 0.16666666666666666, 0.16666666666666666],
        [-0.16666666666666666, 0.16666666666666666, 0.16666666666666666],
        [-0.5, -0.5, 0.16666666666666666],
        [-0.16666666666666666, -0.5, 0.16666666666666666],
        [-0.5, -0.5, -0.5]
    ]
    kweights = [0.07407407407407407,
        0.2222222222222222,
        0.2222222222222222,
        0.2222222222222222,
        0.2222222222222222,
        0.037037037037037035]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)

    size, shift = [3, 3, 2], false
    kcoords = [
        [0.0, 0.0, 0.0],
        [0.3333333333333333, 0.0, 0.0],
        [0.3333333333333333, 0.3333333333333333, 0.0],
        [-0.3333333333333333, 0.3333333333333333, 0.0],
        [0.0, 0.0, -0.5],
        [0.3333333333333333, 0.0, -0.5],
        [0.3333333333333333, 0.3333333333333333, -0.5],
        [-0.3333333333333333, 0.3333333333333333, -0.5],
    ]
    kweights = [
        0.05555555555555555,
        0.2222222222222222,
        0.1111111111111111,
        0.1111111111111111,
        0.05555555555555555,
        0.2222222222222222,
        0.1111111111111111,
        0.1111111111111111,
    ]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)

    size, shift = [3, 3, 2], true
    kcoords = [
        [0.16666666666666666, 0.16666666666666666, 0.25],
        [-0.5, 0.16666666666666666, 0.25],
        [-0.16666666666666666, 0.16666666666666666, 0.25],
        [-0.5, -0.5, 0.25],
        [-0.16666666666666666, -0.5, 0.25],
        [-0.16666666666666666, -0.16666666666666666, 0.25],
    ]
    kweights = [0.1111111111111111,
        0.2222222222222222,
        0.2222222222222222,
        0.1111111111111111,
        0.2222222222222222,
        0.1111111111111111
    ]
    test(3, silicon.lattice, silicon.atoms, silicon.positions, size, shift;
        kcoords=kcoords, kweights=kweights)

    ########### test 2D #############
    lattice = [[1.0 0.0]; [0.0 1.0]]
    atoms = [1,]
    pos = [[0.0, 0.0],]
    size, shift = [2, 2], false
    kcoords = [
        [0.0, 0.0],
        [-0.5, 0.0],
        [-0.5, -0.5]
    ]

    kweights = [] #TODO: not calculated

    meshmap, brmesh = test(2, lattice, atoms, pos, size, shift;
        kcoords=kcoords, kweights=kweights)
    # for i in meshmap.irreducible_indices
    #     println(brmesh[i])
    # end
end