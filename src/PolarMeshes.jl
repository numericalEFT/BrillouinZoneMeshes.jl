
# this file is included in BZMeshes
using ..CoordinateTransformations

export PolarMesh, Angular

const Angular = Union{Polar,Spherical}
const _polar2cart = CartesianFromPolar()
const _spherical2cart = CartesianFromSpherical()
const _cart2polar = PolarFromCartesian()
const _cart2spherical = SphericalFromCartesian()

struct PolarMesh{T,DIM,MT} <: AbstractMesh{T,DIM}
    br::Brillouin{T,DIM}
    mesh::MT # actual mesh. assume order as (r,θ,ϕ...)
end

Base.length(mesh::PolarMesh) = length(mesh.mesh)
Base.size(mesh::PolarMesh) = size(mesh.mesh)
Base.size(mesh::PolarMesh, I::Int) = size(mesh.mesh, I)

# getindex return cartesian results to be consistent with general mesh convention
function Base.getindex(mesh::PolarMesh{T,2,MT}, i::Int, j::Int) where {T,MT}
    return _polar2cart(Polar(getindex(mesh.mesh, i, j)...))
end
function Base.getindex(mesh::PolarMesh{T,3,MT}, i::Int, j::Int, k::Int) where {T,MT}
    return _spherical2cart(Spherical(getindex(mesh.mesh, i, j, k)...))
end
function Base.getindex(mesh::PolarMesh, I::Int)
    return Base.getindex(mesh, _ind2inds(size(mesh), I)...)
end
# provide getindex which return angular results
function Base.getindex(mesh::PolarMesh{T,2,MT}, ::Type{<:Angular}, i::Int, j::Int) where {T,MT}
    return Polar(getindex(mesh.mesh, i, j)...)
end
function Base.getindex(mesh::PolarMesh{T,3,MT}, ::Type{<:Angular}, i::Int, j::Int, k::Int) where {T,MT}
    return Spherical(getindex(mesh.mesh, i, j, k)...)
end
function Base.getindex(mesh::PolarMesh, T::Type{<:Angular}, I::Int)
    return Base.getindex(mesh, T, _ind2inds(size(mesh), I)...)
end
