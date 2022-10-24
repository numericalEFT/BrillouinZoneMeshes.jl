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

# We represent then the BZ as a set of irreducible points `kpoints`,
# and a set of weights `kweights` (summing to 1). The value of
# observables is given by a weighted sum over the irreducible kpoints,
# plus a symmetrization operation (which depends on the particular way
# the observable transforms under the symmetry).

# There is by decreasing cardinality
# - The group of symmetry operations of the lattice
# - The group of symmetry operations of the crystal (model.symmetries)
# - The group of symmetry operations of the crystal that preserves the BZ mesh (basis.symmetries)

# See https://juliamolsim.github.io/DFTK.jl/stable/advanced/symmetries for details.

@doc raw"""
Return the ``k``-point symmetry operations associated to a lattice and atoms.
"""
function symmetry_operations(lattice, atoms, positions, magnetic_moments=[];
    tol_symmetry=SYMMETRY_TOLERANCE)
    @assert length(atoms) == length(positions)
    atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]
    Ws, ws = spglib_get_symmetry(lattice, atom_groups, positions, magnetic_moments; tol_symmetry)
    [SymOp(W, w) for (W, w) in zip(Ws, ws)]
end

# Approximate in; can be performance-critical, so we optimize in case of rationals
is_approx_in_(x::AbstractArray{<:Rational}, X) = any(isequal(x), X)
is_approx_in_(x::AbstractArray{T}, X) where {T} = any(y -> isapprox(x, y; atol=sqrt(eps(T))), X)

"""
Filter out the symmetry operations that don't respect the symmetries of the discrete BZ grid
"""
function symmetries_preserving_kgrid(symmetries, kcoords)
    kcoords_normalized = normalize_kpoint_coordinate.(kcoords)
    function preserves_grid(symop)
        all(is_approx_in_(normalize_kpoint_coordinate(symop.S * k), kcoords_normalized)
            for k in kcoords_normalized)
    end
    filter(preserves_grid, symmetries)
end

"""
Filter out the symmetry operations that don't respect the symmetries of the discrete real-space grid
"""
function symmetries_preserving_rgrid(symmetries, fft_size)
    is_in_grid(r) =
        all(zip(r, fft_size)) do (ri, size)
            abs(ri * size - round(ri * size)) / size ≤ SYMMETRY_TOLERANCE
        end

    onehot3(i) = (x = zeros(Bool, 3); x[i] = true; Vec3(x))
    function preserves_grid(symop)
        all(is_in_grid(symop.W * onehot3(i) .// fft_size[i] + symop.w) for i = 1:3)
    end

    filter(preserves_grid, symmetries)
end

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

""""
Convert a `basis` into one that doesn't use BZ symmetry.
This is mainly useful for debug purposes (e.g. in cases we don't want to
bother thinking about symmetries).
"""
function unfold_bz(basis::PlaneWaveBasis)
    if length(basis.symmetries) == 1
        return basis
    else
        kcoords = unfold_kcoords(basis.kcoords_global, basis.symmetries)
        return PlaneWaveBasis(basis.model,
            basis.Ecut, basis.fft_size, basis.variational,
            kcoords, [1 / length(kcoords) for _ in kcoords],
            basis.kgrid, basis.kshift,
            basis.symmetries_respect_rgrid, basis.comm_kpts)
    end
end

# find where in the irreducible basis `basis_irred` the k-point `kpt_unfolded` is handled
function unfold_mapping(basis_irred, kpt_unfolded)
    for ik_irred = 1:length(basis_irred.kpoints)
        kpt_irred = basis_irred.kpoints[ik_irred]
        for symop in basis_irred.symmetries
            Sk_irred = normalize_kpoint_coordinate(symop.S * kpt_irred.coordinate)
            k_unfolded = normalize_kpoint_coordinate(kpt_unfolded.coordinate)
            if (Sk_irred ≈ k_unfolded) && (kpt_unfolded.spin == kpt_irred.spin)
                return ik_irred, symop
            end
        end
    end
    error("Invalid unfolding of BZ")
end

function unfold_array_(basis_irred, basis_unfolded, data, is_ψ)
    if basis_irred == basis_unfolded
        return data
    end
    if !(basis_irred.comm_kpts == basis_irred.comm_kpts == MPI.COMM_WORLD)
        error("Brillouin zone symmetry unfolding not supported with MPI yet")
    end
    data_unfolded = similar(data, length(basis_unfolded.kpoints))
    for ik_unfolded in 1:length(basis_unfolded.kpoints)
        kpt_unfolded = basis_unfolded.kpoints[ik_unfolded]
        ik_irred, symop = unfold_mapping(basis_irred, kpt_unfolded)
        if is_ψ
            # transform ψ_k from data into ψ_Sk in data_unfolded
            kunfold_coord = kpt_unfolded.coordinate
            @assert normalize_kpoint_coordinate(kunfold_coord) ≈ kunfold_coord
            _, ψSk = apply_symop(symop, basis_irred,
                basis_irred.kpoints[ik_irred], data[ik_irred])
            data_unfolded[ik_unfolded] = ψSk
        else
            # simple copy
            data_unfolded[ik_unfolded] = data[ik_irred]
        end
    end
    data_unfolded
end

function unfold_bz(scfres)
    basis_unfolded = unfold_bz(scfres.basis)
    ψ = unfold_array_(scfres.basis, basis_unfolded, scfres.ψ, true)
    eigenvalues = unfold_array_(scfres.basis, basis_unfolded, scfres.eigenvalues, false)
    occupation = unfold_array_(scfres.basis, basis_unfolded, scfres.occupation, false)
    E, ham = energy_hamiltonian(basis_unfolded, ψ, occupation;
        scfres.ρ, eigenvalues, scfres.εF)
    @assert E.total ≈ scfres.energies.total
    new_scfres = (; basis=basis_unfolded, ψ, ham, eigenvalues, occupation)
    merge(scfres, new_scfres)
end

""""
Convert a `kcoords` into one that doesn't use BZ symmetry.
This is mainly useful for debug purposes (e.g. in cases we don't want to
bother thinking about symmetries).
"""
function unfold_kcoords(kcoords, symmetries)
    all_kcoords = [normalize_kpoint_coordinate(symop.S * kcoord)
                   for kcoord in kcoords, symop in symmetries]

    # the above multiplications introduce an error
    unique(all_kcoords) do k
        digits = ceil(Int, -log10(SYMMETRY_TOLERANCE))
        normalize_kpoint_coordinate(round.(k; digits))
    end
end
