using Brillouin, PlotlyJS
using Brillouin:
    AVec,
    CARTESIAN
import Brillouin:
    basis
using Brillouin.KPaths: Bravais
using StaticArrays
using BrillouinZoneMeshes
using PyCall
const PySpatial = PyNULL()
using BrillouinZoneMeshes.LinearAlgebra
include("default_colors.jl")
include("plotlyjs_wignerseitz.jl")
function wignerseitz_ext(basis::AVec{<:SVector{D,<:Real}};
    merge::Bool=false,
    Nmax::Integer=3, cut=Nmax / 2) where {D}
    # bringing in SciPy's Spatial module (for `Voronoi` and `ConvexHull`)
    if PyCall.conda
        copy!(PySpatial, pyimport_conda("scipy.spatial", "scipy"))
    else
        copy!(PySpatial, pyimport_e("scipy.spatial"))
    end
    ispynull(PySpatial) && @warn("scipy python package not found. " *
                                 "WignerSeitz.wignerseitz is nonfunctional.")

    # "supercell" lattice of G-vectors
    Ns = -Nmax:Nmax
    lattice = Vector{SVector{D,Float64}}(undef, length(Ns)^D)
    idx_cntr = 0
    idx_center = 0
    idx_cartesian = CartesianIndices(ntuple(_ -> Ns, Val(D)))
    for (idx, I) in enumerate(idx_cartesian)
        V = basis[1] * I[1]
        D >= 2 && (V += basis[2] * I[2])
        D == 3 && (V += basis[3] * I[3])
        lattice[idx] = V
        iszero(I) && (idx_center = idx)
    end
    ispynull(PySpatial) && error("You need to install scipy for wignerseitz to work.")
    vor = PySpatial.Voronoi(lattice) # voronoi tesselation of lattice
    clist = Cell{D}[]
    idx_center_final = 0
    # grab all the vertices of the central voronoi region enclosing origo
    # for idx_cntr in 1:length(vor.point_region)
    for idx_cntr in 1:length(vor.point_region)
        if sqrt(sum((idx_cartesian[idx_cntr][i] - idx_cartesian[idx_center][i])^2 for i in 1:D)) < cut
            verts_cntr =  # NB: offsets by 1 due to Julia 1-based vs. Python 0-based indexing
                [vor.vertices[idx+1, :] for idx in vor.regions[vor.point_region[idx_cntr]+1] if !isnothing(idx) && idx > 0]
            # get convex hull of central vertices
            if length(verts_cntr) > D + 1
                hull = PySpatial.ConvexHull(verts_cntr)
                c = WignerSeitz.convert_to_cell(hull, basis)
                c = WignerSeitz.reorient_normals!(c)

                # return either raw simplices or "merged" polygons
                if merge || D > 2
                    latticize!(WignerSeitz.merge_coplanar!(c))
                else
                    latticize!(c)
                end
                push!(clist, c)
                if idx_cntr == idx_center
                    idx_center_final = length(clist)
                end
            end
        end
    end
    return clist, idx_center_final
end

function wignerseitz_ext(basis::AVec{<:AVec{<:Real}}; kwargs...)
    D = length(first(basis))
    all(V -> length(V) == D, @view basis[2:end]) || error(DomainError(basis, "provided `basis` must have identical dimensions"))

    return wignerseitz_ext(SVector{D,Float64}.(basis); kwargs...)
end


# Wigner-Seitz cells visualization
# 2D
#Rs = [[1.0, 0.0], [-0.5, √3/2]]
#N, DIM = 4, 2
#origin = [0.0 0.0]
#lattice =Matrix([2 0; 1 sqrt(3)]')
#msize = (3,3)
#br = BZMeshes.Brillouin(lattice = lattice)

# latvec = [1 0; 0 1]'
# 3D

lattice = [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
atoms = [1, 1]
positions = [ones(3) / 8, -ones(3) / 8]
br = BZMeshes.Brillouin(lattice=lattice, atoms=atoms, positions=positions)
#lattice = [[1.0 0.0 0.0]; [0.0 1.0 0.0]; [0.0 0.0 1.0]]
#br = BZMeshes.Brillouin(lattice=lattice)
msize = (4, 4, 4)

bzmesh = UniformBZMesh(br=br, size=msize)

meshmap = MeshMap(bzmesh)

latvec = mapslices(x -> [x], lattice_vector(bzmesh), dims=1)[:]
clist, idx_center = wignerseitz_ext(latvec, cut=1)

P = plot(clist, idx_center)

fullmesh = [bzmesh[i] for i in 1:length(bzmesh)]
reducedmesh = [bzmesh[i] for i in meshmap.irreducible_indices]

addtraces!(P, scatter3d(x=[r[1] for r in fullmesh], y=[r[2] for r in fullmesh], z=[r[3] for r in fullmesh], mode="markers", marker=attr(size=3)))
addtraces!(P, scatter3d(x=[r[1] for r in reducedmesh], y=[r[2] for r in reducedmesh], z=[r[3] for r in reducedmesh], mode="markers", marker=attr(size=3)))

n = length(P.plot.data) # number of traces, the last two are the full and reduced meshes
allinvis, vis1, vis2 = [true for i in 1:n], [true for i in 1:n], [true for i in 1:n]
allinvis[n-1:n] .= false # [true, ..., true, false, false], hide last two traces
vis1[n] = false # [true, ..., true, true, false], hide the last trace
vis2[n-1] = false # [true, ..., true, false, true], hide the second last trace

d = Dict(
    :updatemenus => [
        attr(
            buttons=[
                attr(
                    args=[
                        # 33, attr(visible=true)
                        "visible", vis1
                    ],
                    label="Complete Mesh",
                    method="restyle"
                ),
                attr(
                    args=[
                        # 33, attr(visible=true)
                        "visible", vis2
                    ],
                    label="Reduced Mesh",
                    method="restyle"
                ),
                attr(
                    args=[
                        # 33, attr(visible=false)
                        "visible", allinvis
                    ],
                    label="No Mesh",
                    method="restyle"
                )
            ],
            direction="down",
            pad_r=10,
            pad_t=10,
            showactive=true,
            x=0.1,
            xanchor="left",
            y=1.1,
            yanchor="top"
        )
    ]
)
relayout!(P, d)

display(P)
readline()

