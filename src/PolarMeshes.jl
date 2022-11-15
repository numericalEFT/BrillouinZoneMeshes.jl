
# this file is included in BZMeshes
include("utilities/coordinatesystems.jl")
using .Coordinates

export PolarMesh, Angular, RescaledGrid, RescaledLogDensedGrid

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

struct PolarMesh{T,DIM,MT<:ProdMesh} <: AbstractMesh{T,DIM}
    cell::Cell{T,DIM}
    mesh::MT # actual mesh. assume order as (r,θ,ϕ...)
    volume::T
end

function PolarMesh(cell::Cell{T,2}, mesh::MT) where {T,MT}
    vol = 0.0
    for j in 1:size(mesh)[2]
        for i in 1:size(mesh)[1]
            r1, r2 = AbstractMeshes.interval(mesh.grids[j], i)
            vol += T(0.5) * (r2^2 - r1^2) * volume(mesh.mesh, j)
        end
    end
    return PolarMesh{T,2,MT}(cell, mesh, vol)
end
function PolarMesh(cell::Cell{T,3}, mesh::MT) where {T,MT}
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
    return PolarMesh{T,3,MT}(cell, mesh, vol)
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

# volume of PolarMesh is different from volume of ProdMesh inside
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

BaseMesh.lattice_vector(mesh::PolarMesh) = mesh.cell.recip_lattice
BaseMesh.inv_lattice_vector(mesh::PolarMesh) = mesh.cell.inv_recip_lattice
BaseMesh.lattice_vector(mesh::PolarMesh, i::Int) = Cells.get_latvec(mesh.cell.recip_lattice, i)
BaseMesh.inv_lattice_vector(mesh::PolarMesh, i::Int) = Cells.get_latvec(mesh.cell.inv_recip_lattice, i)
BaseMesh.cell_volume(mesh::PolarMesh) = mesh.cell.recip_cell_volume

###
# Generate PolarMesh with rescaled log densed grid
###

struct RescaledGrid{T,GT,FT,IFT} <: AbstractGrid{T}
    # rescaled grid with function func()
    bound::SVector{2,T}
    size::Int
    grid::Vector{T} # grid == func.(basegrid.grid)

    #additional info
    basegrid::GT
    func::FT
    invfunc::IFT
    function RescaledGrid(basegrid::GT,
        func::FT, invfunc::IFT) where {GT,FT,IFT}

        T = eltype(GT)
        bound = func.(basegrid.bound)
        size = length(basegrid)
        grid = func.(basegrid.grid)

        return new{T,GT,FT,IFT}(bound, size, grid, basegrid, func, invfunc)
    end
end

# many interface could simply inherit from AbstractGrid
# the following need new implementation
Base.floor(grid::RescaledGrid, x) = floor(grid.basegrid, grid.invfunc(x))
CompositeGrids.Interp.locate(grid::RescaledGrid, x) = CompositeGrids.Interp.locate(grid.basegrid, grid.invfunc(x))

function radial_rescale(; grid::AbstractGrid, DIM::Int)
    if DIM == 2
        func = sqrt
        invfunc = x -> x^2
    elseif DIM == 3
        func = cbrt
        invfunc = x -> x^3
    else
        error("DIM=$DIM not implemented!")
    end

    return RescaledGrid(grid, func, invfunc)
end

function find_kFermi(dispersion, angle; kinit=0.0, krange=nothing)
    if isnothing(krange)
        k0 = kinit
    else
        k0 = krange
    end
    if length(angle) == 1
        return find_zero(k -> dispersion(BZMeshes._polar2cart(Polar(k, angle...))), k0)
    elseif length(angle) == 2
        return find_zero(k -> dispersion(BZMeshes._spherical2cart(Spherical(k, angle...))), k0)
    else
        error("dimension $(length(angle)+1) not implemented!")
    end
end

function RescaledLogDensedGrid(type, bound, densepoints, Nlog, minterval, Nbase, DIM)
    if DIM == 2
        func = sqrt
        invfunc = x -> x^2
    elseif DIM == 3
        func = cbrt
        invfunc = x -> x^3
    else
        error("DIM=$DIM not implemented!")
    end
    rbound = invfunc.(bound)
    rdp = invfunc.(densepoints)
    g = CompositeGrid.LogDensedGrid(type, rbound, rdp, Nlog, minterval, Nbase)
    return radial_rescale(grid=g, DIM=DIM)
end

function kF_densed_kgrids(; dispersion,
    anglemesh,
    bound,
    basegridtype=:cheb,
    Nloggrid=3,
    minterval=0.01,
    Nbasegrid=2,
    DIM=2)
    # assume dispersion==0 has one root for each angle
    k1 = find_kFermi(dispersion, anglemesh[1]; krange=bound)
    # g1 = CompositeGrid.LogDensedGrid(basegridtype, bound, [k1,], Nloggrid, minterval, Nbasegrid)
    g1 = RescaledLogDensedGrid(basegridtype, bound, [k1,], Nloggrid, minterval, Nbasegrid, DIM)
    grids = [g1,]
    for i in 2:length(anglemesh)
        kF = find_kFermi(dispersion, anglemesh[i]; krange=bound)
        # g = CompositeGrid.LogDensedGrid(basegridtype, bound, [kF,], Nloggrid, minterval, Nbasegrid)
        g = RescaledLogDensedGrid(basegridtype, bound, [kF,], Nloggrid, minterval, Nbasegrid, DIM)
        push!(grids, g)
    end
    return grids
end


function BZMeshes.PolarMesh(; dispersion, anglemesh, cell, kmax,
    kwargs...)

    DIM = size(cell.lattice, 1)
    bound = [0.0, kmax]
    println(typeof(dispersion))
    grids = kF_densed_kgrids(; dispersion=dispersion, anglemesh=anglemesh, bound=bound, DIM=DIM, kwargs...)
    cm = ProdMesh(anglemesh, grids)
    pm = PolarMesh(cell, cm)
    return pm
end
