
function MeshMaps.MeshMap(model::UniformKMeshSym.Model, dim; Ecut, kgrid_size::Vector{Int}, kshift=zeros(dim))
    symmetries = model.symmetries

    # kcoords_nosym, kweights_nosym, nosymmetris = UniformKMeshSym.bzmesh_uniform(kgrid, kshift=kshift)
    kcoords, kweights, symmetries = UniformKMeshSym.bzmesh_ir_wedge(kgrid_size, symmetries; kshift=kshift)
    all_kcoords = UniformKMeshSym.unfold_kcoords(kcoords, symmetries)
    # println(kcoords)
    # println(all_kcoords)
    Nk = reduce(*, kgrid_size)
    kmap = zeros(Int, Nk)
    inv_kmap = Dict{Int,Vector{Int}}()
    for kpoint in kcoords
        all_kpoint = UniformKMeshSym.unfold_kcoords([kpoint,], symmetries)
    end

    @assert length(all_kcoords) == reduce(*, kgrid_size)

    return all_kcoords, kcoords
end

