using Unitful
using UnitfulAtomic
using BrillouinZoneMeshes
using DFTK

# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.0]  # Silicon lattice vectors
    [1 0 1.0]  # specified column by column
    [1 1 0.0]];

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
basis = PlaneWaveBasis(model; Ecut, kgrid, kshift=[1 // 2, 1 // 2, 1 // 2])

br = BrillouinZoneMeshes.Model.Brillouin(lattice=austrip.(lattice), atoms=[1, 1], positions=positions)
brmesh = BrillouinZoneMeshes.UniformBZMesh(br=br, size=tuple(kgrid...), shift=1 // 2)

meshmap = BrillouinZoneMeshes.reduced_uniform_meshmap(br, kgrid_size=kgrid, kshift=true)
println(meshmap.irreducible_indices)
for ind in meshmap.irreducible_indices
    # println(AbstractMeshes._ind2inds(tuple(kgrid...), ind))
    println(brmesh[AbstractMeshes._ind2inds(tuple(kgrid...), ind)...])
end