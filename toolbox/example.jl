# using PlotlyJS
# using SymmetryReduceBZ
using BrillouinZoneMeshes
include("Visualization.jl")
# Wigner-Seitz cells visualization
# 2D
#N, DIM = 4, 2
#origin = [0.0 0.0]
#lattice = Matrix([2 0; 1 sqrt(3)]')
#br = BZMeshes.Brillouin(lattice = lattice)

# latvec = [1 0; 0 1]'
# 3D

lattice = [[0 0.5 0.5]; [0.5 0 0.5]; [0.5 0.5 0.0]]
#atoms = [1,1]
#positions = [ones(3) / 8, -ones(3) / 8]

#lattice = [[1.0 0.0 0.0]; [0.0 1.0 0.0]; [0.0 0.0 1.0]]
atoms = [1]
positions = [zeros(3)]
#positions = [zeros(2)]

atom_pos = hvcat(size(positions, 1), positions...)
# ibzformat = "convex hull"
# coordinates = "Cartesian"
# makeprim = true
# convention = "ordinary"

br = BZMeshes.Cell(lattice=lattice, atoms=atoms, positions=positions)
#br = BZMeshes.Brillouin(lattice=lattice)
#print(br.positions, br.lattice, br.atoms)
#msize = (3, 3)
#bzmesh = UniformBZMesh(cell=br, size=msize, shift=[false, false])

msize = (4, 4, 4)
bzmesh = UniformBZMesh(cell=br, size=msize, shift=[true, true, true])
P = Visualization.plotBZ(bzmesh, reduce_to_wigner = false, ex_radius = 1.5,show_reduce = true)

display(P)
readline()


