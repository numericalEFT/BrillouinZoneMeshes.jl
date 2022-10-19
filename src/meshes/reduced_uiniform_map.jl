function _kcoords2ind(kcoord, kgrid_size)
    kidx = [Int(kvec * kgrid_size[d]) + kgrid_size[d] - 1 for (d, kvec) in enumerate(kcoord)]
    klinearidx = BaseMesh._inds2ind(tuple(kgrid_size...), tuple(kidx...))
    return klinearidx
end

function reduced_uniform_meshmap(model::UniformKMeshSym.Model, dim; kgrid_size::Vector{Int})
    symmetries = model.symmetries

    kcoords, kweights, symmetries = UniformKMeshSym.bzmesh_ir_wedge(kgrid_size, symmetries; kshift=zeros(Int, dim))
    all_kcoords = UniformKMeshSym.unfold_kcoords(kcoords, symmetries)
    Nk = reduce(*, kgrid_size)
    @assert length(all_kcoords) == Nk
    kindices = []
    kmap = zeros(Int, Nk)
    inv_kmap = Dict{Int,Vector{Int}}()
    count = 0
    for kpoint in kcoords
        all_kpoint = UniformKMeshSym.unfold_kcoords([kpoint,], symmetries)
        k0ind = _kcoords2ind(kpoint, kgrid_size)
        push!(kindices, k0ind)
        inv_kmap[k0ind] = []
        count += length(all_kpoint)
        for k in all_kpoint
            kind = _kcoords2ind(k, kgrid_size)
            kmap[kind] = k0ind
            push!(inv_kmap[k0ind], kind)
        end
    end

    @assert count == Nk

    return MeshMaps.MeshMap(kindices, kmap, inv_kmap)
end

