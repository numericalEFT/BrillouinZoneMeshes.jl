module Visualization

# this module provides visualization of meshes

using ..BaseMesh
# using ..TreeMeshes
using ..MeshMaps
using Brillouin.WignerSeitz
using Brillouin.KPaths
using PlotlyJS
using Brillouin:
    AVec,
    CARTESIAN,
    cartesianize,
    reduce_to_wignerseitz,
    setting,
    basis
using StaticArrays
#using BrillouinZoneMeshes
#using Brillouin.WignerSeitz.PyCall
using PyCall
const PySpatial = PyNULL()
using BrillouinZoneMeshes.LinearAlgebra
using SymmetryReduceBZ.Symmetry: calc_ibz, inhull, calc_pointgroup, complete_orbit
import SymmetryReduceBZ.Utilities: get_uniquefacets
import QHull
include("default_colors.jl")
include("plotlyjs_wignerseitz.jl")
include("cluster.jl")
# plot mesh as interactive widget via Plotly
export plotBZ

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

function wignerseitz_ext(basis::Matrix; kwargs...)
    latvec = mapslices(x -> [x], basis, dims=1)[:] #Convert matrix into arrays of lattice vectors
    return wignerseitz_ext(latvec; kwargs...)
end

function reduce_to_wignerseitz_ext(v, latvec, bzmesh)
    return cartesianize(reduce_to_wignerseitz(inv_lattice_vector(bzmesh) * v, latvec), latvec)
end


function convert_to_cell(hull::QHull.Chull{Float64}, basis::AVec{<:SVector{D,<:Real}}) where {D}
    # The QHull julia package is used by SymmetryReduceBZ, whereas Brillouin package use Qhull from
    # SciPy instead. Therefore we have to reload this function for julia QHull object.
    vs′ = hull.points         # vertices
    simp = hvcat(size(hull.simplices, 1), hull.simplices...)'
    fs′ = simp # faces

    vs = SVector{D,Float64}.(eachrow(vs′))
    fs = Vector{Int}.(eachrow(fs′))
    return Cell(vs, fs, SVector{D,SVector{D,Float64}}(basis), Ref(CARTESIAN))
end

function convert_to_cell(hull::QHull.Chull{Float64}, basis::Matrix)  
    latvec = mapslices(x -> [x], basis, dims=1)[:] #Convert matrix into arrays of lattice vectors
    return convert_to_cell(hull,  SVector{length(latvec[1]),Float64}.(latvec))
end


function _add_traces!(P,mesh,DIM)
    # P is the plot object
    if DIM == 2
        addtraces!(P, scatter(x=[r[1] for r in mesh], y=[r[2] for r in mesh], mode="markers", marker=attr(size=5)))
    elseif DIM ==3
        addtraces!(P, scatter3d(x=[r[1] for r in mesh], y=[r[2] for r in mesh], z=[r[3] for r in mesh], mode="markers", marker=attr(size=3)))
    else
        error("Dimension must be 2 or 3!")
    end
end

function plotBZ(bzmesh::AbstractUniformMesh{T,DIM}; ex_radius = 1, show_reduce=false, reduce_to_wigner = true,display = nothing , ibzformat = "convex hull", coordinates = "Cartesian",
                makeprim = true, convention = "ordinary"
                ) where {T,DIM}
    # Shows all extended Brillouin Zone within the ex_radius.
    meshmap = MeshMap(bzmesh)
    recip_lattice = lattice_vector(bzmesh)
    clist, idx_center = wignerseitz_ext(recip_lattice, cut=ex_radius)

    fullmesh = []
    for i in 1:length(bzmesh)
        v = bzmesh[i]
        if reduce_to_wigner
            vv = reduce_to_wignerseitz_ext(v, latvec, bzmesh)
            for (gi, g) in enumerate(fullmesh)
                @assert (vv ≈ g) == false "$(bzmesh[i]) and $(bzmesh[gi]) are the same point $(vv) and $(g)"
            end
            push!(fullmesh, vv)
        else
            push!(fullmesh, v)
        end
    end
    # fullmesh = [reduce_to_wignerseitz_ext(bzmesh[i], latvec, bzmesh) for i in 1:length(bzmesh)]

    if show_reduce
        atoms = bzmesh.cell.atoms
        positions = bzmesh.cell.positions
        lattice = bzmesh.cell.lattice
        atom_pos = hvcat(size(positions, 1), positions...)
        #generate irreducible brillouin zone
        ibz = calc_ibz(lattice ./ 2 / π, atoms, atom_pos, coordinates, ibzformat,
                       makeprim, convention)
        #The 1/2π here is due to the convention difference between packages. SymmetryReduceBZ has a \dot a_recip = 1 instead of 2π
        #c = convert_to_cell(ibz, SVector{length(latvec[1]),Float64}.(latvec))
        c = convert_to_cell(ibz,recip_lattice)
        c = WignerSeitz.reorient_normals!(c)
        latticize!(c) # Do not merge coplanar triangles for reduced mesh. This must be done within plot.
        reducedmesh = []
        pointgroup = calc_pointgroup(Matrix{Float64}(lattice_vector(bzmesh)))
        points = [ibz.points[i, :] for i in 1:size(ibz.points, 1)]
        for idx in meshmap.irreducible_indices
            # sym_points = [fullmesh[pidx] for pidx in meshmap.inv_map[idx]]
            orbit = complete_orbit(fullmesh[idx], Matrix{Float64}.(pointgroup))
            orbit = [orbit[:, i] for i in 1:size(orbit, 2)]
            # println(size(orbit))
            point = get_closest_point(points, orbit)
            # println(sym_points, " -> ", point)
            push!(reducedmesh, point)
        end

        P = plot(clist, idx_center, ibz=c)
        #addtraces!(P, scatter(x=[r[1] for r in reducedmesh], y=[r[2] for r in reducedmesh], mode="markers", marker=attr(size=5)))
        _add_traces!(P,reducedmesh, DIM)
    else
        P = plot(clist, idx_center)
    end
    #addtraces!(P, scatter(x=[r[1] for r in fullmesh], y=[r[2] for r in fullmesh], mode="markers", marker=attr(size=5)))
    _add_traces!(P, fullmesh, DIM)

    n = length(P.plot.data) # number of traces, the last two are the full and reduced meshes
    allinvis, vis1, vis2 = [true for i in 1:n], [true for i in 1:n], [true for i in 1:n]
    allinvis[n-1:n] .= false # [true, ..., true, false, false], hide last two traces
    vis1[n] = false # [true, ..., true, true, false], hide the last trace
    vis2[n-1] = false # [true, ..., true, false, true], hide the second last trace

    if isnothing(display)
        display = Dict(
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
    end
    relayout!(P, display)
    return P
end



end
