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
meshmap = MeshMap(bzmesh)

recip_lattice = lattice_vector(bzmesh)
latvec = mapslices(x -> [x], recip_lattice, dims=1)[:]

#generate irreducible brillouin zone
# ibz = calc_ibz(lattice ./ 2 / π, atoms, atom_pos, coordinates, ibzformat,
#     makeprim, convention)
#The 1/2π here is due to the convention difference between packages. SymmetryReduceBZ has a \dot a_recip = 1 instead of 2π

# c = convert_to_cell(ibz, SVector{length(latvec[1]),Float64}.(latvec))
# c = WignerSeitz.reorient_normals!(c)
# latticize!(c) # Do not merge coplanar triangles for reduced mesh. This must be done within plot.


# reducedmesh = [fullmesh[i] for i in meshmap.irreducible_indices]

# reducedmesh = []
# pointgroup = calc_pointgroup(Matrix{Float64}(lattice_vector(bzmesh)))
# points = [ibz.points[i, :] for i in 1:size(ibz.points, 1)]
# for idx in meshmap.irreducible_indices
#     # sym_points = [fullmesh[pidx] for pidx in meshmap.inv_map[idx]]
#     orbit = complete_orbit(fullmesh[idx], Matrix{Float64}.(pointgroup))
#     orbit = [orbit[:, i] for i in 1:size(orbit, 2)]
#     # println(size(orbit))
#     point = get_closest_point(points, orbit)
#     # println(sym_points, " -> ", point)
#     push!(reducedmesh, point)
# end

# remappedmesh = reducedmesh


#addtraces!(P, scatter3d(x=[r[1] for r in fullmesh], y=[r[2] for r in fullmesh], z=[r[3] for r in fullmesh], mode="markers", marker=attr(size=3)))
#addtraces!(P, scatter3d(x=[r[1] for r in reducedmesh], y=[r[2] for r in reducedmesh], z=[r[3] for r in reducedmesh], mode="markers", marker=attr(size=3)))

P = Visualization.plotBZ(bzmesh, reduce_to_wigner = false, ex_radius = 1.5,show_reduce = true)

display(P)
readline()


