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

# tools useful for AbstractMesh
@generated function _inds2ind(size::NTuple{DIM,Int}, I) where {DIM}
    ex = :(I[DIM] - 1)
    for i = (DIM-1):-1:1
        ex = :(I[$i] - 1 + size[$i] * $ex)
    end
    return :($ex + 1)
end

@generated function _ind2inds(size::NTuple{DIM,Int}, I::Int) where {DIM}
    inds, quotient = :((I - 1) % size[1] + 1), :((I - 1) ÷ size[1])
    for i = 2:DIM-1
        inds, quotient = :($inds..., $quotient % size[$i] + 1), :($quotient ÷ size[$i])
    end
    inds = :($inds..., $quotient + 1)
    return :(SVector{DIM,Int}($inds))
end

end