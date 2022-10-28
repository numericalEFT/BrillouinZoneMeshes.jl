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
    meshmap = BrillouinZoneMeshes.uniform_meshmap(brmesh)
    # println(brmehs.)
    for i in meshmap.irreducible_indices
        println(brmesh.mesh.inv_lattice * brmesh[i], ": ", meshmap.inv_map)
        println(brmesh.br.inv_lattice * brmesh[i], ": ", meshmap.inv_map)
        # println(brmesh.br.inv_recip_lattice * brmesh[i], ": ", meshmap.inv_map)
    end
    
    res = 0.0
    for (i, ind) in enumerate(meshmap.irreducible_indices)
        # println(brmesh.mesh.inv_lattice * brmesh[i], ": ", meshmap.inv_map)
        weight = volume(brmesh, ind)*length(meshmap.inv_map[ind])
        res +=  
    end

end