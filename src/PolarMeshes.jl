
# this file is included in BZMeshes
include("utilities/coordinatesystems.jl")
using .Coordinates

export PolarMesh, Angular

const Angular = Union{Polar,Spherical}

# use coordinate systems copied from CoordinateTransformations
# changed some convention thus different from their package
# our ϕ ∈ [-π, π] is the angle from x-axis
# our θ ∈ [-π/2, π/2] is the angle from xy-plane
const _polar2cart = CartesianFromPolar()
const _spherical2cart = CartesianFromSpherical()
const _cart2polar = PolarFromCartesian()
const _cart2spherical = SphericalFromCartesian()

_extract(r::Polar{T,A}) where {T,A} = SVector{2,T}(r.r, r.ϕ)
_extract(r::Spherical{T,A}) where {T,A} = SVector{3,T}(r.r, r.θ, r.ϕ)

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
# call looks like mesh[Angular, i, j] -> Polar(r,θ)
function Base.getindex(mesh::PolarMesh{T,2,MT}, ::Type{<:Angular}, i::Int, j::Int) where {T,MT}
    return Polar(getindex(mesh.mesh, i, j)...)
end
function Base.getindex(mesh::PolarMesh{T,3,MT}, ::Type{<:Angular}, i::Int, j::Int, k::Int) where {T,MT}
    return Spherical(getindex(mesh.mesh, i, j, k)...)
end
function Base.getindex(mesh::PolarMesh, T::Type{<:Angular}, I::Int)
    return Base.getindex(mesh, T, _ind2inds(size(mesh), I)...)
end

function AbstractMeshes.locate(mesh::PolarMesh, r::Angular)
    return AbstractMeshes.locate(mesh.mesh, _extract(r))
end
function AbstractMeshes.locate(mesh::PolarMesh{T,2,MT}, x::AbstractVector) where {T,MT}
    return AbstractMeshes.locate(mesh, _cart2polar(x))
end
function AbstractMeshes.locate(mesh::PolarMesh{T,3,MT}, x::AbstractVector) where {T,MT}
    return AbstractMeshes.locate(mesh, _cart2spherical(x))
end