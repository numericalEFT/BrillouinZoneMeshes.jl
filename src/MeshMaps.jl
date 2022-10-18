module MeshMaps
# symmetry reduce map for meshes

using ..TreeMeshes
using ..BaseMesh
"""
    struct MeshMap

Mapping from full mesh to irreducible mesh. The reduction is mainly available via symmetry operations.

# Parameters:
- `map`: mapping from full mesh to irreducible mesh. When i is the index of a point in the full mesh, map[i] is the corresponding index in reduced mesh
- `reduced_length`: the length of reduced mesh
- `inv_map`: inverse of map. When i is the index of a point in the reduced mesh, inv_map[i] is a list of all points corresponding to this point in the full mesh
"""
struct MeshMap
    irreducible_indices::Vector{Int}
    map::Vector{Int}
    inv_map::Dict{Int,Vector{Int}}

    # function MeshMap(map::Vector{Int})
    #     reduced_length = maximum(map)
    #     inv_map = Vector{Vector{Int}}(undef, reduced_length)
    #     for (i, ind) in enumerate(map)
    #         # i is index of full mesh, ind is index of reduced mesh
    #         if isassigned(inv_map, ind)
    #             push!(inv_map[ind], i)
    #         else
    #             inv_map[ind] = Vector{Int}([i,])
    #         end
    #     end

    #     return new(map, reduced_length, inv_map)
    # end
end

# TODO: constructors that generate map for specific type of mesh and symmetry

## TODO: 1st step: symmetry reduce for M-P mesh(centered uniform mesh)



"""
    struct ReducedMesh{MT}

Map-reduced mesh constructed from mesh::MT with symmetry reduction.

# Parameters:
- `mesh`: bare mesh from which the reduced mesh constructed
- `meshmap`: map from mesh to the reduced mesh
"""
struct ReducedMesh{T,DIM,MT<:AbstractMesh{T,DIM}} <: AbstractMesh{T,DIM}
    mesh::MT
    meshmap::MeshMap
end

# TODO: implement AbstractMesh interface
# including AbstractArray interface and locate/volume functions





######################################################
## LEGACY CODE BELOW
####################################################

export SymMap, MappedData

# TODO: add ReducedMesh: bare mesh + SymMap

function _find_in(x, arr::AbstractArray; atol=1e-6, rtol=1e-6)
    # return index if in, return 0 otherwise
    for yi in 1:length(arr)
        y = arr[yi]
        if isapprox(x, y, atol=atol, rtol=rtol)
            return yi
        end
    end

    return 0
end

# rename to MeshMap
struct SymMap{T,N}
    map::Vector{Int}
    reduced_length::Int
    _vals::Vector{T}
    inv_map::Vector{Vector{Int}}
    # TODO: store symmetry operation

    function SymMap{T}(tg::TreeGrid, density; atol=1e-6, rtol=1e-6) where {T}
        map = zeros(Int, length(tg))
        reduced_vals = []
        inv_map = []
        for pi in 1:length(tg)
            # println(pi, " ", p)
            p = tg[pi]
            val = density(p)
            # println(val)
            pos = _find_in(val, reduced_vals; atol=atol, rtol=rtol)
            if pos == 0
                push!(reduced_vals, val)
                push!(inv_map, [pi,])
                map[pi] = length(reduced_vals)
            else
                push!(inv_map[pos], pi)
                map[pi] = pos
            end
        end

        return new{T,length(tg)}(map, length(reduced_vals), reduced_vals, inv_map)
    end
end

struct MappedData{T,N} <: AbstractArray{T,N}
    smap::SymMap{T,N}
    data::Vector{T}

    function MappedData(smap::SymMap{T,N}) where {T,N}
        data = zeros(T, smap.reduced_length)
        return new{T,N}(smap, data)
    end
end

Base.length(md::MappedData) = length(md.smap.map)
Base.size(md::MappedData) = (length(md),)
# index and iterator
Base.getindex(md::MappedData, i) = md.data[md.smap.map[i]]
function Base.setindex!(md::MappedData, x, i)
    md.data[md.smap.map[i]] = x
end
Base.firstindex(md::MappedData) = 1
Base.lastindex(md::MappedData) = length(tg)

Base.iterate(md::MappedData) = (md[1], 1)
Base.iterate(md::MappedData, state) = (state >= length(md)) ? nothing : (md[state+1], state + 1)

end