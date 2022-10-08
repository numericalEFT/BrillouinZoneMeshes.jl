module PolarMeshes

using ..StaticArrays
using ..CompositeGrids
using ..BaseMesh

export PolarMesh, AngleMesh

_polar2cartesian(r, θ, ϕ) = SVector{3,Float}{r * sin(θ) * cos(ϕ),r * sin(θ) * sin(ϕ),r * cos(θ)}
_polar2cartesian(r, θ) = SVector{2,Float}(r * cos(θ), r * sin(θ))

struct AngleMesh{DIM,AT} <: AbstractMesh{DIM}
    angle_grids::AT # AT is supposed to be a Tuple or Vector of angle grids
    dims::NTuple{DIM,Int}
end

Base.length(mesh::AngleMesh) = prod(mesh.dims)
Base.size(mesh::AngleMesh) = mesh.dims
Base.size(mesh::AngleMesh, I::Int) = mesh.dims[I]

function Base.getindex(mesh::AngleMesh{DIM,AT}, inds...) where {DIM,AT}
end

struct PolarMesh{DIM,RG} <: AbstractMesh{DIM}
    radial_grids::Vector{RG}
    angle_mesh::AngleMesh
end

end