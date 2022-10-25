module MeshMaps
# symmetry reduce map for meshes

using ..AbstractMeshes
using ..Model
using ..TreeMeshes
using ..BaseMesh
using ..BZMeshes

export MeshMap, ReducedBZMesh

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
    inv_indices::Dict{Int,Int}

    function MeshMap(map::Vector{Int})
        irreducible_indices = Vector{Int}([])
        inv_map = Dict{Int,Vector{Int}}([])
        for (i, ind) in enumerate(map)
            if !(ind in irreducible_indices)
                push!(irreducible_indices, ind)
                push!(inv_map, (ind => [i,]))
            else
                push!(inv_map[ind], i)
            end
        end
        inv_indices = Dict{Int,Int}([])
        for (i, irind) in enumerate(irreducible_indices)
            push!(inv_indices, (irind => i))
        end
        return new(irreducible_indices, map, inv_map, inv_indices)
    end
end

# TODO: constructors that generate map for specific type of mesh and symmetry
MeshMap(mesh::AbstractMesh) = error("Map reduce not defined for $(typeof(mesh))!")

Base.length(mm::MeshMap) = length(mm.irreducible_indices)
Base.size(mm::MeshMap) = (length(mm),)
Base.getindex(mm::MeshMap, I::Int) = mm.map[I]

# locate corresponding index in reduced mesh for index of full mesh
AbstractMeshes.locate(mm::MeshMap, I::Int) = mm.inv_indices[mm.map[I]]
_foldnumber(mm::MeshMap, I::Int) = length(mm.inv_map[I])

## TODO: 1st step: symmetry reduce for M-P mesh(centered uniform mesh)


"""
    struct ReducedMesh{MT}

Map-reduced mesh constructed from mesh::MT with symmetry reduction.

# Parameters:
- `mesh`: bare mesh from which the reduced mesh constructed
- `meshmap`: map from mesh to the reduced mesh
"""
struct ReducedBZMesh{T,DIM,MT<:AbstractMesh{T,DIM}} <: AbstractMesh{T,DIM}
    mesh::MT
    meshmap::MeshMap
end

# length and size return REDUCED length of mesh
Base.length(mesh::ReducedBZMesh) = length(mesh.meshmap)
Base.size(mesh::ReducedBZMesh) = size(mesh.meshmap)

# linear index only go through REDUCED points
Base.getindex(mesh::ReducedBZMesh, I::Int) = mesh.mesh[mesh.meshmap.irreducible_indices[I]]
# cartesian index return points of FULL mesh
Base.getindex(mesh::ReducedBZMesh, inds...) = Base.getindex(mesh.mesh, inds...)

AbstractMeshes.locate(mesh::ReducedBZMesh, x) = locate(mesh.meshmap, locate(mesh.mesh, x))
AbstractMeshes.volume(mesh::ReducedBZMesh) = volume(mesh.mesh)
AbstractMeshes.volume(mesh::ReducedBZMesh, I::Int) = volume(mesh.mesh, mesh.meshmap.irreducible_indices[I]) * _foldnumber(mesh.meshmap, mesh.meshmap.irreducible_indices[I])

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