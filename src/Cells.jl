module Cells

using ..PointSymmetry
using ..StaticArrays
using ..LinearAlgebra
import ..showfieldln
using Printf

export Cell, get_latvec

"""
Compute the inverse of the lattice. Require lattice to be square matrix
"""
function _compute_inverse_lattice(lattice::AbstractMatrix{T}) where {T}
    @assert size(lattice, 2) == size(lattice, 2)
    return inv(lattice)
end

"""
Compute the reciprocal lattice.
We use the convention that the reciprocal lattice is the set of G vectors such
that G ⋅ R ∈ 2π ℤ for all R in the lattice.
"""
function _compute_recip_lattice(lattice::AbstractMatrix{T}) where {T}
    return 2T(π) * _compute_inverse_lattice(lattice')
end

"""
Return I-th lattice vector of lattice.
Lattice vectors are specified column-wise in lattice::Matrix.
"""
function get_latvec(lattice::AbstractMatrix{T}, I::Int) where {T}
    return view(lattice, :, I)
end

"""
    struct Brillouin{T, DIM}

Container storing information of Brillouin zone. Including lattice vector, reciprocal lattice vector and their inverse;
volume of unit cell and reciprocal unit cell; G vectors for extended Brillouin zone.

Lattice parameters `lattice` are given by a ``3×3`` matrix with floating point values,
where ``𝐚``, ``𝐛``, and ``𝐜`` are given as __columns__.

# Parameters:
- `T`: type of data
- `DIM`: dimension of the Brillouin zone

# Members:
- `lattice`: lattice vector
- `recip_lattice`: reciprocal lattice vector
- `inv_lattice`: inverse of lattice vector
- `inv_recip_lattice`: inverse of reciprocal lattice vector

- `cell_volume`: volume of lattice unit cell
- `recip_cell_volume`: volume of reciprocal lattice unit cell

- `atoms`: list of integers representing atom types
- `positions`:  # positions[i] is the location of atoms[i] in fract. coords
- `atom_groups``:  atoms[i] == atoms[j] for all i, j in atom_group[α]

- `G_vector`: a list of G vectors in extended Brillouin zone
"""
struct Cell{T,DIM}
    # Lattice and reciprocal lattice vectors in columns
    lattice::Matrix{T}
    recip_lattice::Matrix{T}
    # Useful for conversions between cartesian and reduced coordinates
    inv_lattice::Matrix{T}
    inv_recip_lattice::Matrix{T}

    cell_volume::T
    recip_cell_volume::T

    atoms::Vector{Int}
    positions::Vector{SVector{DIM,T}}  # positions[i] is the location of atoms[i] in fract. coords
    atom_groups::Vector{Vector{Int}}  # atoms[i] == atoms[j] for all i, j in atom_group[α]


    # collection of all allowed G vectors
    G_vector::Vector{SVector{DIM,Int}}
    symmetry::Vector{PointSymmetry.SymOp}
end

function Cell(;
    lattice::Matrix{T},
    atoms::AbstractVector{Int}=Vector{Int}([1,]),
    positions=[zeros(size(lattice, 1)),],
    G_vector=nothing) where {T}

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
    cell_volume = abs(det(lattice))
    recip_cell_volume = abs(det(recip_lattice))

    # G vector default (0,0,0)
    if isnothing(G_vector)
        G_vector = [SVector{DIM,Int}(zeros(DIM)),]
    end

    cell = Cell{T,DIM}(lattice, recip_lattice, inv_lattice, inv_recip_lattice, cell_volume, recip_cell_volume, atoms, positions, atom_groups, G_vector, [])
    for symop in default_symmetries(cell)
        push!(cell.symmetry, symop)
    end
    return cell
end

function get_latvec(br::Cell, I::Int; isrecip=true)
    if isrecip
        return get_latvec(br.recip_lattice, I)
    else
        return get_latvec(br.lattice, I)
    end
end

"""
    function standard_cell(;
        dtype=Float64,
        lattice::AbstractMatrix,
        atoms::AbstractVector=[1,],
        positions::AbstractVector{<:AbstractVector}=[zeros(dtype, size(lattice, 1)),],
        primitive=true,
        correct_symmetry=true,
        # magnetic_moments=[],
        tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE,
        G_vector=nothing)
    
Returns a `Cell` object with crystallographic conventional cell according to the International Table of
Crystallography Vol A (ITA) in case `primitive=false`. If `primitive=true`
the primitive lattice is returned in the convention of the reference work of
Cracknell, Davies, Miller, and Love (CDML). Of note this has minor differences to
the primitive setting choice made in the ITA.
"""
function standard_cell(;
    dtype=Float64,
    lattice::AbstractMatrix,
    atoms::AbstractVector=[1,],
    positions::AbstractVector{<:AbstractVector}=[zeros(dtype, size(lattice, 1)),],
    primitive=true,
    correct_symmetry=true,
    # magnetic_moments=[],
    tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE,
    G_vector=nothing)


    DIM = size(lattice, 1)

    # @assert DIM == 3 "Only 3D lattices are supported."

    # Atoms and terms
    if length(atoms) != length(positions)
        error("Length of atoms and positions vectors need to agree.")
    end
    atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]


    _lattice, _positions = PointSymmetry._make3D(lattice, positions)

    magnetic_moments = []
    cell, _ = PointSymmetry.spglib_cell(_lattice, atom_groups, _positions, magnetic_moments)
    std_cell = PointSymmetry.standardize_cell(cell; to_primitive=primitive, symprec=tol_symmetry,
        no_idealize=!correct_symmetry)

    _lattice = Matrix{dtype}(std_cell.lattice)
    _types = std_cell.types
    _atoms = Int.(_types)
    _positions = Vector{dtype}.(std_cell.positions)
    # magnetic_moments = normalize_magnetic_moment.(std_cell.magmoms)

    ### truncate 3D convention to DIM
    lattice = zeros(dtype, DIM, DIM)
    lattice[1:DIM, 1:DIM] .= _lattice[1:DIM, 1:DIM]
    positions = []
    for i in eachindex(_positions)
        push!(positions, _positions[i][1:DIM])
    end

    ## make sure the truncation doesn't cause error
    __lattice, __positions = PointSymmetry._make3D(lattice, positions)
    @assert __lattice ≈ _lattice "lattice $_lattice is truncated to $__lattice"
    @assert __positions ≈ _positions "position $_positions is truncated to $__position"

    return Cell(;
        lattice=lattice,
        atoms=_atoms,
        positions=positions,
        G_vector=G_vector
    )
end

"""
Default logic to determine the symmetry operations to be used in the model.
"""
function default_symmetries(model::Cell{T,DIM}
    ; tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE) where {T,DIM}

    lattice = zeros(T, 3, 3)
    for i in 1:DIM
        lattice[i, 1:DIM] = model.lattice[i, 1:DIM]
    end
    for i in DIM+1:3
        lattice[i, i] = 1
    end
    positions = [zeros(T, 3) for i in 1:length(model.atoms)]
    for ai in eachindex(model.atoms)
        positions[ai][1:DIM] = model.positions[ai][1:DIM]
    end
    atoms = model.atoms

    # Standard case from here on:
    if length(positions) != length(atoms)
        error("Length of atoms and positions vectors need to agree.")
    end

    return PointSymmetry.symmetry_operations(lattice, atoms, positions; tol_symmetry)
end


# TODO: Add the following helper functions
# #=
# There are two types of quantities, depending on how they transform under change of coordinates.

# Positions transform with the lattice: r_cart = lattice * r_red. We term them vectors.

# Linear forms on vectors (anything that appears in an expression f⋅r) transform
# with the inverse lattice transpose: if f_cart ⋅ r_cart = f_red ⋅ r_red, then
# f_cart = lattice' \ f_red. We term them covectors.
# Examples of covectors are forces.

# Reciprocal vectors are a special case: they are covectors, but conventionally have an
# additional factor of 2π in their definition, so they transform rather with 2π times the
# inverse lattice transpose: q_cart = 2π lattice' \ q_red = recip_lattice * q_red.
# =#

# TODO: define following funcs for cells

vector_frac_to_cart(cell::Cell, rfrac) = cell.lattice * rfrac
vector_cart_to_frac(cell::Cell, rcart) = cell.inv_lattice * rcart
covector_frac_to_cart(cell::Cell, ffrac) = cell.inv_lattice' * ffrac
covector_cart_to_frac(cell::Cell, fcart) = cell.lattice' * fcart
recip_vector_frac_to_cart(cell::Cell, qfrac) = cell.recip_lattice * qfrac
recip_vector_cart_to_frac(cell::Cell, qcart) = cell.inv_recip_lattice * qcart

export vector_cart_to_frac, vector_frac_to_cart
export covector_cart_to_frac, covector_frac_to_cart
export recip_vector_cart_to_frac, recip_vector_frac_to_cart
# #=
# Transformations on vectors and covectors are matrices and comatrices.

# Consider two covectors f and g related by a transformation matrix B. In reduced
# coordinates g_red = B_red f_red and in cartesian coordinates we want g_cart = B_cart f_cart.
# From g_cart = L⁻ᵀ g_red = L⁻ᵀ B_red f_red = L⁻ᵀ B_red Lᵀ f_cart, we see B_cart = L⁻ᵀ B_red Lᵀ.

# Similarly for two vectors r and s with s_red = A_red r_red and s_cart = A_cart r_cart:
# s_cart = L s_red = L A_red r_red = L A_red L⁻¹ r_cart, thus A_cart = L A_red L⁻¹.

# Examples of matrices are the symmetries in real space (W)
# Examples of comatrices are the symmetries in reciprocal space (S)
# =#
# matrix_red_to_cart(model::Brillouin, Ared) = model.lattice * Ared * model.inv_lattice
# matrix_cart_to_red(model::Brillouin, Acart) = model.inv_lattice * Acart * model.lattice
# comatrix_red_to_cart(model::Brillouin, Bred) = model.inv_lattice' * Bred * model.lattice'
# comatrix_cart_to_red(model::Brillouin, Bcart) = model.lattice' * Bcart * model.inv_lattice'

# Implementation of the show function for Cell

function Base.show(io::IO, model::Cell{T,DIM}) where {T,DIM}
    print(io, "Cell (", DIM, "D) with lattice vectors $(model.lattice))")
end

function Base.show(io::IO, ::MIME"text/plain", model::Cell{T,DIM}) where {T,DIM}
    println(io, "Cell (", DIM, "D)")
    for i = 1:DIM
        header = i == 1 ? "lattice" : ""
        if DIM == 1
            showfieldln(io, header, (@sprintf "[%-10.6g, ]" model.lattice[i, :]...))
        elseif DIM == 2
            showfieldln(io, header, (@sprintf "[%-10.6g, %-10.6g]" model.lattice[i, :]...))
        elseif DIM == 3
            showfieldln(io, header, (@sprintf "[%-10.6g, %-10.6g, %-10.6g]" model.lattice[i, :]...))
        else
            showfieldln(io, header, ("$(model.lattice[i, :]...)"))
        end
    end
    showfieldln(io, "unit cell volume", @sprintf "%.5g" model.cell_volume)

    if !isempty(model.atoms)
        println(io)
        showfieldln(io, "atoms", model.atoms)
        for (i, el) in enumerate(model.atoms)
            header = "$i-th atom $el position"
            showfieldln(io, header, model.positions[i])
        end
        # showfieldln(io, "positions", model.positions)
    end
end

end
