using Unitful
using UnitfulAtomic
using DFTK
using StaticArrays

@testset "DFTK mesh interface" begin

    function test_grid(kgrid, kshift::Bool)
        println("Test grid: ", kgrid, " with a shift ", kshift)
        _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]
        # 1. Define lattice and atomic positions
        a = 5.431u"angstrom"          # Silicon lattice constant
        lattice = a / 2 * [[0 1 1.0]  # Silicon lattice vectors
            [1 0 1.0]  # specified column by column
            [1 1 0.0]]

        positions = [ones(3) / 8, -ones(3) / 8]

        ######## DFTK interface ########################
        # 2. Select model and basis
        Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))

        # Specify type and positions of atoms
        atoms = [Si, Si]
        model = model_LDA(lattice, atoms, positions)

        kgrid = [4, 4, 4]     # k-point grid (Regular Monkhorst-Pack grid)
        Ecut = 7              # kinetic energy cutoff
        # Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
        basis = PlaneWaveBasis(model; Ecut, kgrid, kshift=_kshift)

        br = BrillouinZoneMeshes.Model.Brillouin(lattice=austrip.(lattice), atoms=[1, 1], positions=positions)
        brmesh = BrillouinZoneMeshes.UniformBZMesh(br=br, size=tuple(kgrid...), shift=_kshift[1])

        meshmap = BrillouinZoneMeshes.reduced_uniform_meshmap(br, kgrid_size=kgrid, kshift=kshift)
        println(meshmap.irreducible_indices)
        for (i, ind) in enumerate(meshmap.irreducible_indices)
            # println(AbstractMeshes._ind2inds(tuple(kgrid...), ind))
            # println(brmesh[AbstractMeshes._ind2inds(tuple(kgrid...), ind)...])
            kvec = brmesh[AbstractMeshes._ind2inds(tuple(kgrid...), ind)...]
            _kvec = DFTK.recip_vector_red_to_cart(model, basis.kcoords_global[i])
            # println(kvec, " vs ", _kvec)
            @test kvec ≈ _kvec
        end
    end
    test_grid([4, 4, 4], false)
    test_grid([4, 4, 4], true)
    test_grid([4, 4, 5], false)
    test_grid([4, 4, 5], true)

end