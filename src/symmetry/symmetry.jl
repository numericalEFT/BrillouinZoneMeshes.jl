# This file contains functions to handle the symetries.
# The type SymOp is defined in Symop.jl

# A symmetry (W, w) (or (S, τ)) induces a symmetry in the Brillouin zone
# that the Hamiltonian at S k is unitary equivalent to that at k, which we exploit to
# reduce computations. The relationship is
# S = W'
# τ = -W^-1 w
# (valid both in reduced and cartesian coordinates). In our notation
# the rotation matrix W and translation w are such that, for each atom of
# type A at position a, W a + w is also an atom of type A.

# The full (reducible) Brillouin zone is implicitly represented by
# a set of (irreducible) kpoints (see explanation in docs). Each
# irreducible k-point k comes with a list of symmetry operations
# (S, τ) (containing at least the trivial operation (I, 0)), where S
# is a unitary matrix (/!\ in cartesian but not in reduced coordinates)
# and τ a translation vector. The k-point is then used to represent
# implicitly the information at all the kpoints Sk. The relationship
# between the Hamiltonians is
# H_{Sk} = U H_k U*, with
# (Uu)(x) = u(W x + w)
# or in Fourier space
# (Uu)(G) = e^{-i G τ} u(S^-1 G)
# In particular, we can choose the eigenvectors at Sk as u_{Sk} = U u_k

# There is by decreasing cardinality
# - The group of symmetry operations of the lattice
# - The group of symmetry operations of the crystal (cell.symmetry)
# - The group of symmetry operations of the crystal that preserves the BZ mesh (mesh.symmetry)

# See https://juliamolsim.github.io/DFTK.jl/stable/advanced/symmetries for details.

"""Bring ``k``-point coordinates into the range [-0.5, 0.5)"""
function normalize_kpoint_coordinate(x::Real)
    x = x - round(Int, x, RoundNearestTiesUp)
    @assert -0.5 ≤ x < 0.5
    x
end
normalize_kpoint_coordinate(k::AbstractVector) = normalize_kpoint_coordinate.(k)

@doc raw"""
Return the ``k``-point symmetry operations associated to a lattice and atoms.

# Arguments
- lattice :: 3x3 AbstractMatrix: lattice vectors in columns.
- atoms :: Vector of atoms
- positions :: Vector of positions of the atoms in fractional coordinates.
"""
function symmetry_operations(lattice, atoms, positions, magnetic_moments=[];
    tol_symmetry=SYMMETRY_TOLERANCE)
    @assert length(atoms) == length(positions)
    atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]
    Ws, ws = spglib_get_symmetry(lattice, atom_groups, positions, magnetic_moments; tol_symmetry)
    return [SymOp(W, w) for (W, w) in zip(Ws, ws)]
end

# Approximate in; can be performance-critical, so we optimize in case of rationals
is_approx_in_(x::AbstractArray{<:Rational}, X) = any(isequal(x), X)
is_approx_in_(x::AbstractArray{T}, X) where {T} = any(y -> isapprox(x, y; atol=sqrt(eps(T))), X)

# """
# Filter out the symmetry operations that don't respect the symmetries of the discrete BZ grid
# """
# function symmetries_preserving_kgrid(symmetries, kcoords)
#     kcoords_normalized = normalize_kpoint_coordinate.(kcoords)
#     function preserves_grid(symop)
#         all(is_approx_in_(normalize_kpoint_coordinate(symop.S * k), kcoords_normalized)
#             for k in kcoords_normalized)
#     end
#     filter(preserves_grid, symmetries)
# end

@doc raw"""
Apply various standardisations to a lattice and a list of atoms. It uses spglib to detect
symmetries (within `tol_symmetry`), then cleans up the lattice according to the symmetries
(unless `correct_symmetry` is `false`) and returns the resulting standard lattice
and atoms. If `primitive` is `true` (default) the primitive unit cell is returned, else
the conventional unit cell is returned.
"""
function standardize_atoms(lattice, atoms, positions, magnetic_moments=[]; kwargs...)
    @assert length(atoms) == length(positions)
    @assert isempty(magnetic_moments) || (length(atoms) == length(magnetic_moments))
    atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]
    ret = spglib_standardize_cell(lattice, atom_groups, positions, magnetic_moments; kwargs...)
    (; ret.lattice, atoms, ret.positions, ret.magnetic_moments)
end

# """"
# Convert a `kcoords` into one that doesn't use BZ symmetry.
# This is mainly useful for debug purposes (e.g. in cases we don't want to
# bother thinking about symmetries).
# """
# function unfold_kcoords(kcoords, symmetries)
#     all_kcoords = [normalize_kpoint_coordinate(symop.S * kcoord)
#                    for kcoord in kcoords, symop in symmetries]

#     # the above multiplications introduce an error
#     unique(all_kcoords) do k
#         digits = ceil(Int, -log10(SYMMETRY_TOLERANCE))
#         normalize_kpoint_coordinate(round.(k; digits))
#     end
# end

"""
Filter out the kcoords that are irreducible under the given symmetries
"""
function irreducible_kcoord(symmetries::AbstractVector{SymOp}, kcoords::AbstractVector)
    kcoords = _makeVec3.(kcoords)
    kidxmap = [i for i in 1:length(kcoords)]
    for (ki, k) in enumerate(kcoords)
        if kidxmap[ki] != ki # already mapped to a symmetry equivalent
            continue
        end
        symk = [symop.S * k for symop in symmetries]
        for kj in ki+1:length(kcoords)
            if kidxmap[kj] != kj # already mapped to a symmetry equivalent
                continue
            end
            for s in 1:length(symmetries)
                # If the difference between kred and W' * k == W^{-1} * k
                # is only integer in fractional reciprocal-space coordinates, then
                # kred and S' * k are equivalent k-points
                if all(isinteger, kcoords[kj] - symk[s])
                    kidxmap[kj] = ki
                end
            end
        end
    end
    return kidxmap
end

