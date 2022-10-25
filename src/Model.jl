module Model

using ..PointSymmetry
using ..StaticArrays
using ..LinearAlgebra
import ..showfieldln
using Printf

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
    unit_cell_volume = abs(det(lattice))
    recip_cell_volume = abs(det(recip_lattice))

    # G vector default (0,0,0)
    if isnothing(G_vector)
        G_vector = [SVector{DIM,Int}(zeros(DIM)),]
    end

    return Brillouin{T,DIM}(lattice, recip_lattice, inv_lattice, inv_recip_lattice, unit_cell_volume, recip_cell_volume, atoms, positions, atom_groups, G_vector)
end

# normalize_magnetic_moment(::Nothing)::Vec3{Float64} = (0, 0, 0)
# normalize_magnetic_moment(mm::Number)::Vec3{Float64} = (0, 0, mm)
# normalize_magnetic_moment(mm::AbstractVector)::Vec3{Float64} = mm

"""
    function standard_brillouin(;
        dtype=Float64,
        lattice::AbstractMatrix,
        atoms::AbstractVector=[1,],
        positions::AbstractVector{<:AbstractVector}=[zeros(dtype, size(lattice, 1)),],
        primitive=true,
        correct_symmetry=true,
        # magnetic_moments=[],
        tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE,
        G_vector=nothing)
    
Returns a `Brillouin` object with crystallographic conventional cell according to the International Table of
Crystallography Vol A (ITA) in case `primitive=false`. If `primitive=true`
the primitive lattice is returned in the convention of the reference work of
Cracknell, Davies, Miller, and Love (CDML). Of note this has minor differences to
the primitive setting choice made in the ITA.
"""
function standard_brillouin(;
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
    T = dtype

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
    # println(std_cell.positions)
    _positions = Vector{dtype}.(std_cell.positions)
    # magnetic_moments = normalize_magnetic_moment.(std_cell.magmoms)
    # println(_lattice)
    for i in 1:DIM
        lattice[1:DIM, 1:DIM] .= _lattice[1:DIM, 1:DIM]
    end
    for i in eachindex(_positions)
        positions[i][1:DIM] .= _positions[i][1:DIM]
    end
    __lattice, __positions = PointSymmetry._make3D(lattice, positions)
    @assert __lattice ≈ _lattice
    @assert __positions ≈ _positions

    return Brillouin(;
        lattice=lattice,
        atoms=atoms,
        positions=positions,
        G_vector=G_vector)
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
# vector_red_to_cart(model::Brillouin, rred) = model.lattice * rred
# vector_cart_to_red(model::Brillouin, rcart) = model.inv_lattice * rcart
# covector_red_to_cart(model::Brillouin, fred) = model.inv_lattice' * fred
# covector_cart_to_red(model::Brillouin, fcart) = model.lattice' * fcart
# recip_vector_red_to_cart(model::Brillouin, qred) = model.recip_lattice * qred
# recip_vector_cart_to_red(model::Brillouin, qcart) = model.inv_recip_lattice * qcart

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

# Implementation of the show function for Model

function Base.show(io::IO, model::Brillouin{T,DIM}) where {T,DIM}
    print(io, "Brillouin Zone (", DIM, "D) with lattice vectors $(model.lattice))")
end

function Base.show(io::IO, ::MIME"text/plain", model::Brillouin{T,DIM}) where {T,DIM}
    println(io, "Brillouin Zone (", DIM, "D)")
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
    showfieldln(io, "unit cell volume", @sprintf "%.5g" model.unit_cell_volume)

    if !isempty(model.atoms)
        println(io)
        showfieldln(io, "atoms", model.atoms)
        for (i, el) in enumerate(model.atoms)
            header = "atom $el position"
            showfieldln(io, header, model.positions[i])
        end
        # showfieldln(io, "positions", model.positions)
    end
end

end