
# this file is included in BZMeshes
using ..CoordinateTransformations

export PolarMesh, AngleMesh

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

function Base.getindex(mesh::PolarMesh{T,2,MT}, inds...) where {T,MT}
    return _polar2cart(Polar(getindex(mesh.mesh, inds...)))
end
function Base.getindex(mesh::PolarMesh{T,3,MT}, inds...) where {T,MT}
    return _spherical2cart(Polar(getindex(mesh.mesh, inds...)))
end
