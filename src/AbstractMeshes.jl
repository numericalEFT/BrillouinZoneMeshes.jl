module AbstractMeshes

# define interface of AbstractMesh

using ..StaticArrays

export AbstractMesh, locate, volume

# the return value of AbstractMesh should be a SVector{T,DIM}
abstract type AbstractMesh{T,DIM} <: AbstractArray{SVector{T,DIM},DIM} end

Base.IteratorSize(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = Base.HasLength()
Base.IteratorEltype(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = Base.HasEltype()
Base.eltype(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = eltype(T)

# these iterators relies on the interfaces to be implemented below
# and might be implemented in other ways
Base.firstindex(mesh::AbstractMesh) = 1
Base.lastindex(mesh::AbstractMesh) = length(mesh)
Base.iterate(mesh::AbstractMesh) = (mesh[1], 1)
Base.iterate(mesh::AbstractMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

# below are interfaces that should be implemented by concrete types
Base.length(mesh::AbstractMesh) = error("not implemented!")
Base.size(mesh::AbstractMesh) = error("not implemented!")
Base.size(mesh::AbstractMesh, I) = error("not implemented!")

Base.show(io::IO, mesh::AbstractMesh) = error("not implemented!")

Base.getindex(mesh::AbstractMesh, inds...) = error("not implemented!")
Base.getindex(mesh::AbstractMesh, I) = error("not implemented!")

locate(mesh::AbstractMesh, x) = error("not implemented!")
volume(mesh::AbstractMesh) = error("not implemented!")
volume(mesh::AbstractMesh, I) = error("not implemented!")

end