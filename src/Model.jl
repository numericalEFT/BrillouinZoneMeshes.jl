module Model

using ..StaticArrays
using ..LinearAlgebra

export Brillouin

"""
Compute the inverse of the lattice. Require lattice to be square matrix
"""
function _compute_inverse_lattice(lattice::Matrix{T}) where {T}
    @assert size(lattice, 2) == size(lattice, 2)
    return inv(lattice)
end

"""
Compute the reciprocal lattice.
We use the convention that the reciprocal lattice is the set of G vectors such
that G ⋅ R ∈ 2π ℤ for all R in the lattice.
"""
function _compute_recip_lattice(lattice::Matrix{T}) where {T}
    return 2T(π) * _compute_inverse_lattice(lattice)
end

"""
    struct Brillouin{T, DIM}

Container storing information of Brillouin zone. Including lattice vector, reciprocal lattice vector and their inverse;
volume of unit cell and reciprocal unit cell; G vectors for extended Brillouin zone.

# Parameters:
- `T`: type of data
- `DIM`: dimension of the Brillouin zone

# Members:
- `lattice`: lattice vector
- `recip_lattice`: reciprocal lattice vector
- `inv_lattice`: inverse of lattice vector
- `inv_recip_lattice`: inverse of reciprocal lattice vector
- `unit_cell_volume`: volume of lattice unit cell
- `recip_cell_volume`: volume of reciprocal lattice unit cell
- `G_vector`: a list of G vectors in extended Brillouin zone
"""
struct Brillouin{T,DIM}
    # Lattice and reciprocal lattice vectors in columns
    lattice::Matrix{T}
    recip_lattice::Matrix{T}
    # Useful for conversions between cartesian and reduced coordinates
    inv_lattice::Matrix{T}
    inv_recip_lattice::Matrix{T}

    unit_cell_volume::T
    recip_cell_volume::T

    # Particle types (elements) and particle positions and in fractional coordinates.
    # Possibly empty. It's up to the `term_types` to make use of this (or not).
    # `atom_groups` contains the groups of indices into atoms and positions, which
    # point to identical atoms. It is computed automatically on Model construction and may
    # be used to optimise the term instantiation.
    atoms::Vector{Int}
    positions::Vector{SVector{DIM,T}}  # positions[i] is the location of atoms[i] in fract. coords
    atom_groups::Vector{Vector{Int}}  # atoms[i] == atoms[j] for all i, j in atom_group[α]


    # collection of all allowed G vectors
    G_vector::Vector{SVector{DIM,Int}}
end


function Brillouin(;
    lattice::Matrix{T},
    atoms::AbstractVector{Int}=Vector{Int}([]),
    positions=nothing,
    G_vector::Union{Nothing,AbstractVector}=nothing) where {T}

    DIM = size(lattice, 1)

    if isnothing(positions)
        positions = Vector{SVector{DIM,T}}[]
    end

    # Atoms and terms
    if length(atoms) != length(positions)
        error("Length of atoms and positions vectors need to agree.")
    end
    atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]

    # Lattice Vectors
    @assert size(lattice, 1) == size(lattice, 2) "Lattice vector should be given in square matrix!"
    recip_lattice = _compute_recip_lattice(lattice)
    inv_lattice = _compute_inverse_lattice(lattice)
    inv_recip_lattice = _compute_inverse_lattice(recip_lattice)
    unit_cell_volume = abs(det(lattice))
    recip_cell_volume = abs(det(recip_lattice))

    # G vector default (0,0,0)
    if isnothing(G_vector)
        G_vector = [SVector{DIM,Int}(zeros(DIM)),]
    end

    return Brillouin{T,DIM}(lattice, recip_lattice, inv_lattice, inv_recip_lattice, unit_cell_volume, recip_cell_volume, atoms, positions, atom_groups, G_vector)
end

end