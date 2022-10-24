using Unitful
using UnitfulAtomic
using DFTK
using StaticArrays

@testset "DFTK mesh interface" begin

    function test_grid_3D(kgrid, kshift::Bool, dim=3)
        println("Test grid: ", kgrid, " with a shift ", kshift)

        ksize = tuple(kgrid[1:dim]...)
        _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]
        # kshift should be converted between our(VASP) and DFTK convention
        _kshift_dftk = [(iseven(ksize[i]) ? _kshift[i] : _kshift[i] + 1 // 2) for i in 1:length(_kshift)]
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

        # kgrid = [4, 4, 4]     # k-point grid (Regular Monkhorst-Pack grid)
        Ecut = 7              # kinetic energy cutoff
        # Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
        basis = PlaneWaveBasis(model; Ecut, kgrid, kshift=_kshift_dftk)

        br = BrillouinZoneMeshes.Model.Brillouin(lattice=austrip.(lattice)[1:dim, 1:dim], atoms=[1, 1], positions=[positions[1][1:dim], positions[2][1:dim]])
        brmesh = BrillouinZoneMeshes.UniformBZMesh(br=br, size=tuple(kgrid[1:dim]...), shift=_kshift[1])

        meshmap = BrillouinZoneMeshes._reduced_uniform_meshmap(br, kgrid_size=kgrid[1:dim], kshift=kshift)
        println(meshmap.irreducible_indices)
        for (i, ind) in enumerate(meshmap.irreducible_indices)
            # println(AbstractMeshes._ind2inds(tuple(kgrid...), ind))
            # println(brmesh[AbstractMeshes._ind2inds(tuple(kgrid...), ind)...])
            kvec = brmesh[AbstractMeshes._ind2inds(tuple(kgrid[1:dim]...), ind)...]
            _kvec = DFTK.recip_vector_red_to_cart(model, basis.kcoords_global[i])
            # println(kvec, " vs ", _kvec)
            @test kvec ≈ _kvec
        end
    end
    test_grid_3D([4, 4, 4], false, 3)
    test_grid_3D([4, 4, 4], true, 3)
    test_grid_3D([4, 4, 5], false, 3)
    test_grid_3D([4, 4, 5], true, 3)

    function test_grid_2D(kgrid, kshift::Bool, dim=2)
        println("Test grid: ", kgrid, " with a shift ", kshift)
        _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]
        # 1. Define lattice and atomic positions
        a = 5.431u"angstrom"          # Silicon lattice constant
        lattice = a / 2 * [[0 1]  # Silicon lattice vectors
            [1 0]  # specified column by column
        ]

        positions = [zeros(dim),]

        ######## DFTK interface ########################
        # 2. Select model and basis
        Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))

        # Specify type and positions of atoms
        atoms = [Si,]
        model = model_LDA(lattice, atoms, positions)

        # kgrid = [4, 4, 4]     # k-point grid (Regular Monkhorst-Pack grid)
        Ecut = 7              # kinetic energy cutoff
        # Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
        basis = PlaneWaveBasis(model; Ecut, kgrid, kshift=_kshift)

        br = BrillouinZoneMeshes.Model.Brillouin(lattice=austrip.(lattice)[1:dim, 1:dim], atoms=[1, 1], positions=[positions[1][1:dim], positions[2][1:dim]])
        brmesh = BrillouinZoneMeshes.UniformBZMesh(br=br, size=tuple(kgrid[1:dim]...), shift=_kshift[1])

        meshmap = BrillouinZoneMeshes._reduced_uniform_meshmap(br, kgrid_size=kgrid[1:dim], kshift=kshift)
        println(meshmap.irreducible_indices)
        for (i, ind) in enumerate(meshmap.irreducible_indices)
            # println(AbstractMeshes._ind2inds(tuple(kgrid...), ind))
            # println(brmesh[AbstractMeshes._ind2inds(tuple(kgrid...), ind)...])
            kvec = brmesh[AbstractMeshes._ind2inds(tuple(kgrid[1:dim]...), ind)...]
            _kvec = DFTK.recip_vector_red_to_cart(model, basis.kcoords_global[i])
            # println(kvec, " vs ", _kvec)
            @test kvec ≈ _kvec
        end
    end

    test_grid([4, 4, 1], false, 2)
    test_grid([4, 4, 1], true, 2)
    test_grid([4, 5, 1], false, 2)
    test_grid([4, 5, 1], true, 2)

end