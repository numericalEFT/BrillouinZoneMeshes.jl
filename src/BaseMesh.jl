module BaseMesh

using ..StaticArrays
using ..LinearAlgebra
using ..CompositeGrids

using ..BaryCheb
using ..AbstractMeshes
using ..Cells
using ..Cells: get_latvec

export UMesh, ChebMesh
export UniformMesh, BaryChebMesh, CenteredMesh, EdgedMesh

############## Abstract Uniform Mesh #################
"""
    abstract type AbstractUniformMesh{T,DIM} <: AbstractMesh{T,DIM}

Parent type of all uniform meshes. 

Mesh points of a uniform mesh is assumed to be uniformly distributed 
in an parrallelogram area.
Uniform meshes support fractional coordinates.

All concrete implementations of this abstract type are assumed to have
the following fields:
# Required Fields:
- `origin`: the origin(bottom-left point) of the area
- `shift`: the fractional coordinate shift of the mesh. This is useful to reproduce M-P mesh commonly used in DFT.

and have the following methods implemented:
# Required Methods:
- `lattice_vector`: return lattice vector of the represented area
- `inv_lattice_vector`: return inverse of lattice vector
- `cell_volume`: return cell volume
"""
abstract type AbstractUniformMesh{T,DIM} <: AbstractMesh{T,DIM} end
export AbstractUniformMesh
export inv_lattice_vector, lattice_vector, cell_volume
export fractional_coordinates, cartesian_coordinates

function Base.getindex(mesh::AbstractUniformMesh{T,DIM}, inds...) where {T,DIM}
    n = SVector{DIM,Int}(inds)
    return mesh.origin + lattice_vector(mesh) * ((n .- 1 .+ mesh.shift) ./ mesh.size)
end

function Base.getindex(mesh::AbstractUniformMesh, I::Int)
    return Base.getindex(mesh, AbstractMeshes._ind2inds(mesh.size, I)...)
end

function Base.getindex(mesh::AbstractUniformMesh{T,DIM}, ::Type{<:FracCoords}, I::Int) where {T,DIM}
    n = SVector{DIM,Int}(AbstractMeshes._ind2inds(mesh.size, I))
    return inv_lattice_vector(mesh) * mesh.origin + (n .- 1 .+ mesh.shift) ./ mesh.size
end

"""
    function AbstractMeshes.locate(mesh::AbstractUniformMesh{T,DIM}, x) where {T,DIM}

locate mesh point in mesh that is nearest to x. Useful for Monte-Carlo algorithm.
Could also be used for zeroth order interpolation.

# Parameters
- `mesh`: aimed mesh
- `x`: cartesian pos to locate
"""
function AbstractMeshes.locate(mesh::AbstractUniformMesh{T,DIM}, x) where {T,DIM}
    # find index of nearest grid point to the point
    svx = SVector{DIM,T}(x)
    inds = cart_to_frac(mesh, svx - mesh.origin) .* mesh.size .+ 1.5 .- mesh.shift .+ 2 .* eps.(T.(mesh.size))
    indexall = 1
    factor = 1
    indexall += (cycling_floor(inds[1], mesh.size[1]) - 1) * factor
    for i in 2:DIM
        factor *= mesh.size[i-1]
        indexall += (cycling_floor(inds[i], mesh.size[i]) - 1) * factor
    end

    return indexall
end

"""
    function AbstractMeshes.volume(mesh::AbstractUniformMesh, i)

volume represented by mesh point i. When i is omitted return volume of the whole mesh. 
For M-P mesh it's always volume(mesh)/length(mesh), but for others things are more complecated.
Here we assume periodic boundary condition so for all case it's the same.

# Parameters:
- `mesh`: mesh
- `i`: index of mesh point, if ommited return volume of whole mesh
"""
AbstractMeshes.volume(mesh::AbstractUniformMesh) = cell_volume(mesh)
AbstractMeshes.volume(mesh::AbstractUniformMesh, i) = cell_volume(mesh) / length(mesh)

AbstractMeshes.interp(data, mesh::AbstractUniformMesh, x) = data[locate(mesh, x)]

"""
    function AbstractMeshes.integrate(data, mesh::AbstractUniformMesh)

Default integration for uniform meshes. Use zeroth-order integration, 
i.e. average value times volume.

# Parameters:
- `data`: data
- `mesh`: mesh
"""
function AbstractMeshes.integrate(data, mesh::AbstractUniformMesh)
    result = 0.0
    for i in 1:length(mesh)
        result += data[i] * volume(mesh, i)
    end
    return result
end

struct HasLattice <: AbstractMeshes.LatticeStyle end # has mesh.lattice and mesh.inv_lattice
AbstractMeshes.lattice_vector(::HasLattice, mesh) = mesh.lattice
AbstractMeshes.inv_lattice_vector(::HasLattice, mesh) = mesh.inv_lattice
AbstractMeshes.lattice_vector(::HasLattice, mesh, i::Int) = get_latvec(mesh.lattice, i)
AbstractMeshes.inv_lattice_vector(::HasLattice, mesh, i::Int) = get_latvec(mesh.inv_lattice, i)
AbstractMeshes.cell_volume(::HasLattice, mesh) = mesh.cell_volume

############## general purposed uniform mesh #################
"""
    struct UMesh{T,DIM} <: AbstractUniformMesh{T,DIM}

Simplest uniform mesh with lattice/inv_lattice/cell_volume stored.

# Fields:
- `lattice`: lattice vector
- `inv_lattice`: inverse lattice vector
- `cell_volume`: volume of the area represented
- `origin`: the origin(bottom-left point) of the area
- `size`: size of the mesh
- `shift`: the fractional coordinate shift of the mesh. This is useful to reproduce M-P mesh commonly used in DFT.
"""
struct UMesh{T,DIM} <: AbstractUniformMesh{T,DIM}
    lattice::Matrix{T}
    inv_lattice::Matrix{T}
    cell_volume::T

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    shift::SVector{DIM,Rational}
end

"""
    UMesh(;br::Cell{T,DIM}, origin, size, shift)

Construct a UMesh from a brillouin zone with given origin, size and shift.

# Parameters:
- `br`: brillouin zone containing information of the area represented
- `origin`: the origin(bottom-left point) of the area
- `size`: size of the mesh
- `shift`: the fractional coordinate shift of the mesh. This is useful to reproduce M-P mesh commonly used in DFT.
"""
UMesh(;
    br::Cell{T,DIM},
    origin,
    size,
    shift) where {T,DIM} = UMesh{T,DIM}(
    br.recip_lattice,
    br.inv_recip_lattice,
    br.recip_cell_volume,
    SVector{DIM,T}(br.recip_lattice * origin),
    size,
    SVector{DIM,Rational}(shift)
)

AbstractMeshes.LatticeStyle(::Type{<:UMesh{T,DIM}}) where {T,DIM} = HasLattice()

function Base.show(io::IO, mesh::UMesh)
    println("UMesh with $(length(mesh)) mesh points")
end

function cycling_floor(I, N)
    ifloor = (floor(Int, I) + N) % N
    if ifloor == 0
        return N
    else
        return ifloor
    end
end

include("OrthogonalMeshes.jl")

include("ChebMeshes.jl")

end

