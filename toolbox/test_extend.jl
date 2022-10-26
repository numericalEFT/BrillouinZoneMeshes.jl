using Brillouin, PlotlyJS
using Brillouin:
    AVec,
    CARTESIAN
import Brillouin:
    basis
using Brillouin.KPaths: Bravais
using StaticArrays

using PyCall
const PySpatial = PyNULL()

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
# 3D
Rs = [[1.0, 0.0], [-0.5, √3/2]]
#Rs = 2π.*[[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]]
clist,idx_center = wignerseitz_ext(Rs;merge = false)
#P0=plot(clist[idx_center])
#display(P0)
#readline()
v=[1.0,1.0]
vc = cartesianize(v,Rs)
#println(vc[1:1],vc[2:2])
linspace = range(0,1,length=30)
lin2 = collect(Iterators.product(linspace,linspace))
P=plot(clist,idx_center)
addtraces!(P, scatter(x=[vc[1].*lin2[i][1] for i in 1:length(lin2)],y=[vc[2].*lin2[i][2] for i in 1:length(lin2)],mode="markers"))
display(P)
readline()

