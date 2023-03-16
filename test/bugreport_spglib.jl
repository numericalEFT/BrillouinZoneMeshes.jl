# this file fails to run due to some bugs
# a test should be added when the bug is fixed

using BrillouinZoneMeshes
using Test

@testset "spglib" begin
    @testset "N=12" begin
        kgrid = [12, 12, 12]
        lattice = [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]]
        atoms = [1,]
        positions = [[0.0, 0.0, 0.0],]
        cell = Cells.Cell(; lattice=lattice, atoms=atoms, positions=positions)
        bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid, shift=[false, false, false])
        meshmap = MeshMaps.MeshMap(bzmesh)
    end
    @testset "N=20" begin
        kgrid = [20, 20, 20]
        lattice = [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]]
        atoms = [1,]
        positions = [[0.0, 0.0, 0.0],]
        cell = Cells.Cell(; lattice=lattice, atoms=atoms, positions=positions)
        bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid, shift=[false, false, false])
        meshmap = MeshMaps.MeshMap(bzmesh)
    end
end
