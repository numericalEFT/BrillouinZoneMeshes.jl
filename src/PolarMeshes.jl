
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

struct PolarMesh{T,DIM,MT<:CompositeMesh} <: AbstractMesh{T,DIM}
    br::Brillouin{T,DIM}
    mesh::MT # actual mesh. assume order as (r,θ,ϕ...)
    volume::T
end

function PolarMesh(br::Brillouin{T,2}, mesh::MT) where {T,MT}
    vol = 0.0
    for j in 1:size(mesh)[2]
        for i in 1:size(mesh)[1]
            r1, r2 = AbstractMeshes.interval(mesh.grids[j], i)
            vol += T(0.5) * (r2^2 - r1^2) * volume(mesh.mesh, j)
        end
    end
    return PolarMesh{T,2,MT}(br, mesh, vol)
end
function PolarMesh(br::Brillouin{T,3}, mesh::MT) where {T,MT}
    vol = 0.0
    for k in 1:size(mesh)[3]
        for j in 1:size(mesh)[2]
            for i in 1:size(mesh)[1]
                J = Base._sub2ind(size(mesh)[2:3], j, k)
                r1, r2 = AbstractMeshes.interval(mesh.grids[J], i)
                θ1, θ2 = AbstractMeshes.interval(mesh.mesh.grids[k], j)
                # notice that θ ∈ [-π/2,π/2], so integrand is r^2drd(sin(θ))dϕ
                vol += (r2^3 - r1^3) / 3 * (sin(θ2) - sin(θ1)) * volume(mesh.mesh.mesh, k)
            end
        end
    end
    return PolarMesh{T,3,MT}(br, mesh, vol)
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

# volume of PolarMesh is different from volume of CompositeMesh inside
function AbstractMeshes.volume(mesh::PolarMesh)
    return mesh.volume
end
function AbstractMeshes.volume(mesh::PolarMesh{T,2,MT}, I::Int) where {T,MT}
    i, j = AbstractMeshes._ind2inds(size(mesh), I)
    r1, r2 = AbstractMeshes.interval(mesh.mesh.grids[j], i)
    return T(0.5) * (r2^2 - r1^2) * volume(mesh.mesh.mesh, j)
end
function AbstractMeshes.volume(mesh::PolarMesh{T,3,MT}, I::Int) where {T,MT}
    i, j, k = AbstractMeshes._ind2inds(size(mesh), I)
    J = Base._sub2ind(size(mesh)[2:3], j, k)
    r1, r2 = AbstractMeshes.interval(mesh.mesh.grids[J], i)
    θ1, θ2 = AbstractMeshes.interval(mesh.mesh.mesh.grids[k], j)
    return (r2^3 - r1^3) / 3 * (sin(θ2) - sin(θ1)) * volume(mesh.mesh.mesh.mesh, k)
end

BaseMesh.lattice_vector(mesh::PolarMesh) = mesh.br.recip_lattice
BaseMesh.inv_lattice_vector(mesh::PolarMesh) = mesh.br.inv_recip_lattice
BaseMesh.lattice_vector(mesh::PolarMesh, i::Int) = Model.get_latvec(mesh.br.recip_lattice, i)
BaseMesh.inv_lattice_vector(mesh::PolarMesh, i::Int) = Model.get_latvec(mesh.br.inv_recip_lattice, i)
BaseMesh.cell_volume(mesh::PolarMesh) = mesh.br.recip_cell_volume
