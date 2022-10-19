using Unitful
using UnitfulAtomic
using BrillouinZoneMeshes

# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.0]  # Silicon lattice vectors
    [1 0 1.0]  # specified column by column
    [1 1 0.0]];

# Specify type and positions of atoms
atoms = [1, 1]
positions = [ones(3) / 8, -ones(3) / 8]

# 2. Select model and basis
model = BrillouinZoneMeshes.UniformKMeshSym.Model(austrip.(lattice), atoms, positions)
kgrid = [4, 4, 4]     # k-point grid (Regular Monkhorst-Pack grid)
Ecut = 7              # kinetic energy cutoff
# Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
basis = BrillouinZoneMeshes.UniformKMeshSym.PlaneWaveBasis(model; Ecut, kgrid)

meshmap = BrillouinZoneMeshes.reduced_uniform_meshmap(model, 3, kgrid_size=[4, 4, 4])
println(meshmap.irreducible_indices)