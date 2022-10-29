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
            merge::Bool = true,
            Nmax::Integer = 3, cut = Nmax-1) where D
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
    idx_cartesian = CartesianIndices(ntuple(_->Ns, Val(D)))
    for (idx, I) in enumerate(idx_cartesian)
                   V =  basis[1]*I[1]
        D >= 2 && (V += basis[2]*I[2])
        D == 3 && (V += basis[3]*I[3])
        lattice[idx] = V
        iszero(I) && (idx_center = idx)
    end
    ispynull(PySpatial) && error("You need to install scipy for wignerseitz to work.")
    vor = PySpatial.Voronoi(lattice) # voronoi tesselation of lattice
    clist = Cell{D}[]
    idx_center_final =0
    # grab all the vertices of the central voronoi region enclosing origo
    # for idx_cntr in 1:length(vor.point_region)
    for idx_cntr in 1:length(vor.point_region)
        if sqrt(sum((idx_cartesian[idx_cntr][i]-idx_cartesian[idx_center][i])^2 for i in 1:D))<cut
            verts_cntr =  # NB: offsets by 1 due to Julia 1-based vs. Python 0-based indexing
                [vor.vertices[idx+1,:] for idx in vor.regions[vor.point_region[idx_cntr]+1] if !isnothing(idx)&&idx>0]
            # get convex hull of central vertices
            if length(verts_cntr)>D+1
                hull = PySpatial.ConvexHull(verts_cntr)
                c    = WignerSeitz.convert_to_cell(hull, basis)
                c    = WignerSeitz.reorient_normals!(c)

                # return either raw simplices or "merged" polygons
                if merge
                    latticize!(merge_coplanar!(c))
                else
                    latticize!(c)
                end
                push!(clist,c)
                if idx_cntr==idx_center
                    idx_center_final =length(clist)
                end
            end
        end
    end
    return clist, idx_center_final
end

function wignerseitz_ext(basis::AVec{<:AVec{<:Real}}; kwargs...)
    D = length(first(basis))
    all(V->length(V) == D, @view basis[2:end]) || error(DomainError(basis, "provided `basis` must have identical dimensions"))

    return wignerseitz_ext(SVector{D,Float64}.(basis); kwargs...)
end


# Wigner-Seitz cells visualization
# 2D
Rs = [[1.0, 0.0], [-0.5, √3/2]]
N, DIM = 4, 2
origin = [0.0 0.0]
lattice =Matrix([2 0; 1 sqrt(3)]')
msize = (3,3)
# latvec = [1 0; 0 1]'
br = BZMeshes.Brillouin(lattice = lattice)
bzmesh = UniformBZMesh(br=br, size=msize)
latvec = mapslices(x->[x],lattice_vector(bzmesh),dims=1)[:]
clist,idx_center = wignerseitz_ext(latvec;merge = false)
P=plot(clist,idx_center)
addtraces!(P, scatter(x=[bzmesh[i][1] for i in 1:length(bzmesh)],y=[bzmesh[i][2] for i in 1:length(bzmesh)],mode="markers"))
display(P)
readline()

