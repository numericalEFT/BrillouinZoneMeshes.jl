"""
Default logic to determine the symmetry operations to be used in the model.
"""
function default_symmetries(model::Model.Brillouin{T,DIM}
    ; tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE) where {T,DIM}

    lattice = zeros(T, DIM, DIM)
    for i in 1:DIM
        lattice[i, 1:DIM] = model.lattice[i, 1:DIM]
    end
    positions = [zeros(T, DIM) for i in 1:length(model.atoms)]
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

function _kcoords2ind(kcoord, kgrid_size, kshift)
    kidx = [Int(kvec * kgrid_size[d] - kshift[d]) + kgrid_size[d] - 1 for (d, kvec) in enumerate(kcoord)]
    klinearidx = AbstractMeshes._inds2ind(tuple(kgrid_size...), tuple(kidx...))
    return klinearidx
end

# WARNINING: Do not support Gamma_centered: origin=0  !!!
# Monkhorst-Pack: origin=-1/2, consistent with VASP
# to be consistent with DFTK: 
#  - N is even, VASP is the same as DFTK: shift=0 will include Gamma point, shift=1/2 will not
#  - N is odd, VASP is different as DFTK: shift=0 will not include Gamma point, shift=1/2 will
function _reduced_uniform_meshmap(model::Model.Brillouin{T,DIM}, symmetry::Bool=true;
    kgrid_size::Vector{Int}, kshift::Bool=false,
    tol_symmetry=PointSymmetry.SYMMETRY_TOLERANCE
) where {T,DIM}
    # Determine symmetry operations to use
    if symmetry
        symmetries = default_symmetries(model, tol_symmetry=tol_symmetry)
    else
        symmetries = [one(PointSymmetry.SymOp)]
    end
    @assert !isempty(symmetries)  # Identity has to be always present.
    _kgrid_size = ones(Int, DIM)
    _kgrid_size[1:DIM] = kgrid_size[1:DIM]
    _kshift = kshift ? [1 // 2, 1 // 2, 1 // 2] : [0, 0, 0]

    kcoords, kweights, symmetries = PointSymmetry.bzmesh_ir_wedge(_kgrid_size, symmetries; kshift=_kshift)
    all_kcoords = PointSymmetry.unfold_kcoords(kcoords, symmetries)
    Nk = reduce(*, kgrid_size)
    @assert length(all_kcoords) == Nk
    kindices = []
    kmap = zeros(Int, Nk)
    inv_kmap = Dict{Int,Vector{Int}}()
    count = 0
    for kpoint in kcoords
        all_kpoint = PointSymmetry.unfold_kcoords([kpoint,], symmetries)
        k0ind = _kcoords2ind(kpoint, kgrid_size, _kshift)
        push!(kindices, k0ind)
        inv_kmap[k0ind] = []
        count += length(all_kpoint)
        for k in all_kpoint
            kind = _kcoords2ind(k, kgrid_size, _kshift)
            kmap[kind] = k0ind
            push!(inv_kmap[k0ind], kind)
        end
    end

    @assert count == Nk

    return MeshMaps.MeshMap(kindices, kmap, inv_kmap)
end











