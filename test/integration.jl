@testset "Square Lattice" begin

    gap = 0.2
    # x, y, z is from [-π/2, π/2]
    f(x, y, z) = 1 / (1 - cos(x) * cos(y) * cos(z) + gap) / π^3

    ksize = [8, 8, 8]
    lattice = [1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0]
    br = Brillouin(lattice=lattice, atoms=[1,], positions=[[0.0, 0.0, 0.0],])
    brmesh = BZMeshes.Monkhorst_Pack(br=br, size=tuple(ksize...), shift=[0, 0, 0])
    meshmap = MeshMap(brmesh)

    res = 0.0
    for (i, ind) in enumerate(meshmap.irreducible_indices)
        println(inv_lattice_vector(brmesh) * brmesh[i], ": ", length(meshmap.inv_map[ind]))
        weight = AbstractMeshes.volume(brmesh, ind) * length(meshmap.inv_map[ind])
        x, y, z = brmesh[ind]
        res += f(x, y, z) * weight
    end
    println(res)

end