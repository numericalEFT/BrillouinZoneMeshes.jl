using BrillouinZoneMeshes.PointSymmetry: spglib_spacegroup_number, spglib_standardize_cell
using LinearAlgebra
using Test

@testset "spglib" begin
    a = 10.3
    Si = 1
    Ge = 2

    # silicon
    lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
    atoms = [Si, Si]
    positions = [ones(3) / 8, -ones(3) / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 227
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # silicon with Cartesian x coordinates flipped
    lattice = a / 2 * [[0 -1 -1.0]; [1 0 1.0]; [1 1 0.0]]
    atoms = [Si, Si]
    positions = [ones(3) / 8, -ones(3) / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 227
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # silicon with different lattice vectors
    lattice = a / 2 * [[0 1 1.0]; [-1 0 1.0]; [-1 1 0.0]]
    atoms = [Si, Si]
    positions = [[-1, 1, 1] / 8, -[-1, 1, 1] / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 227
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # Zincblende structure
    lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
    atoms = [Si, Ge]
    positions = [ones(3) / 8, -ones(3) / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 216
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # Zincblende structure with Cartesian x coordinates flipped
    lattice = a / 2 * [[0 -1 -1.0]; [1 0 1.0]; [1 1 0.0]]
    atoms = [Si, Ge]
    positions = [ones(3) / 8, -ones(3) / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 216
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # Zincblende structure with different lattice vectors
    lattice = a / 2 * [[0 1 1.0]; [-1 0 1.0]; [-1 1 0.0]]
    atoms = [Si, Ge]
    positions = [[-1, 1, 1] / 8, -[-1, 1, 1] / 8]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 216
    @test spglib_standardize_cell(model).lattice ≈ a * I(3)

    # Hexagonal close packed lattice
    lattice = a * [[1.0 -1 / 2 0.0]; [0.0 sqrt(3) / 2 0.0]; [0.0 0.0 sqrt(8 / 3)]]
    atoms = [Si, Si]
    positions = [[0, 0, 0], [1 / 3, 2 / 3, 1 / 2]]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 194
    @test spglib_standardize_cell(model).lattice ≈ lattice

    # Hexagonal close packed lattice with x coordinates flipped
    lattice = a * [[-1.0 1 / 2 0.0]; [0.0 sqrt(3) / 2 0.0]; [0.0 0.0 sqrt(8 / 3)]]
    atoms = [Si, Si]
    positions = [[0, 0, 0], [1 / 3, 2 / 3, 1 / 2]]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 194
    @test !(spglib_standardize_cell(model).lattice ≈ lattice)

    # Hexagonal close packed lattice with different lattice vectors
    lattice = a * [[-1.0 -1 / 2 0.0]; [0.0 sqrt(3) / 2 0.0]; [0.0 0.0 sqrt(8 / 3)]]
    atoms = [Si, Si]
    positions = [[0, 0, 0], [-1 / 3, 2 / 3, 1 / 2]]
    model = BrillouinZoneMeshes.BaseMesh.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
    @test spglib_spacegroup_number(model) == 194
    @test (spglib_standardize_cell(model).lattice
           ≈
           a * [[1.0 -1 / 2 0.0]; [0.0 sqrt(3) / 2 0.0]; [0.0 0.0 sqrt(8 / 3)]])
end

@testset "Spglib grid index test" begin
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