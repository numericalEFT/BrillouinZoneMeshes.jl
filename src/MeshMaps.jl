module MeshMaps
# symmetry reduce map for meshes

using ..AbstractMeshes
using ..Cells
# using ..TreeMeshes
using ..BaseMesh
import ..showfieldln
# using ..BZMeshes

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

function Base.show(io::IO, mesh::ReducedBZMesh)
    print(io, "Reduced BZ Mesh")
    print(io, ", derived from ", typeof(mesh.mesh))
    print(io, ", reduction = ", length(mesh), "/", length(mesh.mesh))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::ReducedBZMesh)
    println(io, "Reduced BZ Mesh")
    showfieldln(io, "original mesh = ", mesh.mesh)
    showfieldln(io, "map = ", mesh.meshmap)

    println(io)
end

end