
# this file is included in BZMeshes

export PolarMesh, AngleMesh

struct AngleMesh{DIM,AT} <: AbstractMesh{Float64,DIM}
    angle_grids::AT # AT is supposed to be a Tuple or Vector of angle grids
    dims::NTuple{DIM,Int}
end

Base.length(mesh::AngleMesh) = prod(mesh.dims)
Base.size(mesh::AngleMesh) = mesh.dims
Base.size(mesh::AngleMesh, I::Int) = mesh.dims[I]

function Base.getindex(mesh::AngleMesh{DIM,AT}, inds...) where {DIM,AT}
end

struct PolarMesh{DIM,RG} <: AbstractMesh{Float64,DIM}
    radial_grids::Vector{RG}
    angle_mesh::AngleMesh
end

