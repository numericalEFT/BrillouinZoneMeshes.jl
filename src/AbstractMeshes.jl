module AbstractMeshes

# define interface of AbstractMesh

using ..StaticArrays
using ..CompositeGrids

export AbstractCoords, CartCoords, FracCoords, AngularCoords
export dimension
export AbstractMesh, locate, volume#, fractional_coordinates
export FracCoords, frac_to_cart, cart_to_frac
export LatticeStyle, lattice_vector, inv_lattice_vector, cell_volume
export interp, integrate

# the default return value of AbstractMesh should be a SVector{T,DIM}
abstract type AbstractMesh{T,DIM} <: AbstractArray{SVector{T,DIM},DIM} end

# domain
abstract type LatticeStyle end
struct UnknownLattice <: LatticeStyle end # default unknown
abstract type HasCell <: LatticeStyle end # has mesh.cell::Cell
struct BrillouinLattice <: HasCell end # has mesh.cell and use recip lattice
struct BravaisLattice <: HasCell end # has mesh.cell and use lattice

# coordinate types
abstract type AbstractCoords end
struct CartCoords <: AbstractCoords end # default cartesian, stored in SVector
struct FracCoords <: AbstractCoords end # fractional, also in SVector. lattice * frac = cart

Base.IteratorSize(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = Base.HasLength()
Base.IteratorEltype(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = Base.HasEltype()
Base.eltype(::Type{AbstractMesh{T,DIM}}) where {T,DIM} = eltype(T)

# these iterators relies on the interfaces to be implemented below
# and might be implemented in other ways
Base.firstindex(mesh::AbstractMesh) = 1
Base.lastindex(mesh::AbstractMesh) = length(mesh)
Base.iterate(mesh::AbstractMesh) = (mesh[1], 1)
Base.iterate(mesh::AbstractMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

# assume mesh.size::NTuple{DIM,Int}, need implementation otherwise
Base.length(mesh::AbstractMesh) = prod(mesh.size)
Base.size(mesh::AbstractMesh) = mesh.size
Base.size(mesh::AbstractMesh, I) = mesh.size[I]
dimension(mesh::AbstractMesh{T,DIM}) where {T,DIM} = DIM

# below are interfaces that should be implemented by concrete types
Base.show(io::IO, mesh::AbstractMesh) = error("not implemented!")

# getindex dispatch with respect to traits "Coords"
# default should be implemented to return cartesian coords in SVector
# FracCoords and AngularCoords can be given as flags to require corresponding coords
# FracCoords are also in SVector
# while AngularCoords return concrete type Polar<:AngularCoords and Spherical<:AngularCoords

# default, return cartesian coords, should be implemented for all meshes
Base.getindex(mesh::AbstractMesh, inds...) = error("not implemented!")
Base.getindex(mesh::AbstractMesh, I) = error("not implemented!")

# FracCoords, return fractional coordinates in SVector, only apply for meshes with lattice info
Base.getindex(mesh::AbstractMesh, ::Type{<:FracCoords}, inds...) = error("not implemented")
Base.getindex(mesh::AbstractMesh, ::Type{<:FracCoords}, I) = error("not implemented")

# AngularCoords are implemented in PolarMeshes.jl

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

# optional functions

# wrapper of external functions from CompositeGrids
locate(grid::AbstractGrid, x) = CompositeGrids.Interp.locate(grid, x[1])
volume(grid::AbstractGrid) = CompositeGrids.Interp.volume(grid)
volume(grid::AbstractGrid, I) = CompositeGrids.Interp.volume(grid, I)
function interval(grid::AbstractGrid, i::Int)
    # bounds of the inverval around I needed for volume of PolarMesh
    if i != 1 && i != length(grid)
        return (grid[i-1] + grid[i]) / 2, (grid[i+1] + grid[i]) / 2
    elseif i == 1
        return grid.bound[1], (grid[i+1] + grid[i]) / 2
    else
        return (grid[i] + grid[i-1]) / 2, grid.bound[2]
    end
end

# lattice vector information
# abstract type LatticeStyle end # dispatch type for lattice vector functions
# struct UnknownLattice <: LatticeStyle end # no info by default

LatticeStyle(::Type) = UnknownLattice()

lattice_vector(mesh::MT) where {MT} = lattice_vector(LatticeStyle(MT), mesh)
lattice_vector(mesh::MT, i::Int) where {MT} = lattice_vector(LatticeStyle(MT), mesh, i)
inv_lattice_vector(mesh::MT) where {MT} = inv_lattice_vector(LatticeStyle(MT), mesh)
inv_lattice_vector(mesh::MT, i::Int) where {MT} = inv_lattice_vector(LatticeStyle(MT), mesh, i)
cell_volume(mesh::MT) where {MT} = cell_volume(LatticeStyle(MT), mesh)

# by default return error
lattice_vector(::UnknownLattice, mesh) = error("no lattice information!")
lattice_vector(::UnknownLattice, mesh, i::Int) = error("no lattice information!")
inv_lattice_vector(::UnknownLattice, mesh) = error("no lattice information!")
inv_lattice_vector(::UnknownLattice, mesh, i::Int) = error("no lattice information!")
cell_volume(::UnknownLattice, mesh) = error("no lattice information!")

# conversion between fractional and cartesian 
frac_to_cart(mesh::AbstractMesh, frac) = lattice_vector(mesh) * frac
cart_to_frac(mesh::AbstractMesh, cart) = inv_lattice_vector(mesh) * cart

# optional: interp and integrate
integrate(data, mesh::AbstractMesh) = error("not implemented!")
interp(data, mesh::AbstractMesh, x) = error("not implemented!")

# optional: interval. useful for converting from cartesian to angular
interval(mesh::AbstractMesh, I::Int) = error("not implemented!")

end