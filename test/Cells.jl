@testset "Cells" begin
    using BrillouinZoneMeshes.PointSymmetry: spglib_spacegroup_number, spglib_standardize_cell

    @testset "Brillouin" begin
        Brillouin = BrillouinZoneMeshes.Cells.Cell
        get_latvec = BrillouinZoneMeshes.Cells.get_latvec
        # square lattice
        DIM = 2
        lattice = Matrix([1.0 0; 0 1]')
        br = Brillouin(lattice=lattice)
        @test get_latvec(br, 1; isrecip=false) ≈ [1.0, 0]
        @test get_latvec(br, 1) ≈ [2π, 0]
        @test br.inv_lattice .* 2π ≈ br.recip_lattice'
        @test br.cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM

        # helper functions
        x = [1.0, 0.5]
        @test BrillouinZoneMeshes.Cells.vector_frac_to_cart(br, x) == br.lattice * x
        @test BrillouinZoneMeshes.Cells.vector_cart_to_frac(br, x) == br.inv_lattice * x
        @test BrillouinZoneMeshes.Cells.covector_frac_to_cart(br, x) == br.inv_lattice' * x
        @test BrillouinZoneMeshes.Cells.covector_cart_to_frac(br, x) == br.lattice' * x
        @test BrillouinZoneMeshes.Cells.recip_vector_frac_to_cart(br, x) == br.recip_lattice * x
        @test BrillouinZoneMeshes.Cells.recip_vector_cart_to_frac(br, x) == br.inv_recip_lattice * x


        # triagular lattice
        DIM = 2
        lattice = Matrix([2.0 0; 1 sqrt(3)]')
        br = Brillouin(lattice=lattice)
        @test br.inv_lattice .* 2π ≈ br.recip_lattice'
        @test br.cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM
        for i in 1:DIM
            @test dot(get_latvec(br.recip_lattice, i), get_latvec(br.lattice, i)) ≈ 2π
        end
        # 3d testing lattice
        DIM = 3
        lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
        br = Brillouin(lattice=lattice)
        @test br.inv_lattice .* 2π ≈ br.recip_lattice'
        @test br.cell_volume ≈ abs(det(lattice))
        @test br.recip_cell_volume ≈ 1 / abs(det(lattice)) * (2π)^DIM
        for i in 1:DIM
            @test dot(get_latvec(br.recip_lattice, i), get_latvec(br.lattice, i)) ≈ 2π
        end

    end

    @testset "Standard Brillouin" begin
        a = 10.3
        Si = 1
        Ge = 2

        # silicon with Cartesian x coordinates flipped
        lattice = a / 2 * [[0 -1 -1.0]; [1 0 1.0]; [1 1 0.0]]
        atoms = [Si, Si]
        positions = [ones(3) / 8, -ones(3) / 8]
        model = BrillouinZoneMeshes.Cells.standard_cell(lattice=lattice, atoms=atoms, positions=positions, primitive=false)
        println(model)
        display(model)

        @test spglib_spacegroup_number(model) == 227
        @test model.lattice ≈ a * I(3)

        # Zincblende structure with different lattice vectors
        lattice = a / 2 * [[0 1 1.0]; [-1 0 1.0]; [-1 1 0.0]]
        atoms = [Si, Ge]
        positions = [[-1, 1, 1] / 8, -[-1, 1, 1] / 8]
        model = BrillouinZoneMeshes.Cells.standard_cell(lattice=lattice, atoms=atoms, positions=positions, primitive=false)
        @test spglib_spacegroup_number(model) == 216
        @test model.lattice ≈ a * I(3)

        # Two-dimensional example
        lattice = [[1.0 0.0]; [0.0 1.0]]
        atoms = [Si,]
        positions = [[0.0, 0.0],]
        model = BrillouinZoneMeshes.Cells.standard_cell(lattice=lattice, atoms=atoms, positions=positions, primitive=true)
        @test model.lattice ≈ lattice
    end
end
