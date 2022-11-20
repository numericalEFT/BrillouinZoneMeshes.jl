module CompositeMeshes

using ..StaticArrays
using ..LinearAlgebra

using ..AbstractMeshes
using ..AbstractMeshes: _inds2ind, _ind2inds
using ..BaseMesh
using ..CompositeGrids

export CompositeMesh

# CompositeMesh is defined as a panel mesh and a set of submeshes

struct CompositeMesh{T,PM,SM} <: AbstractMesh{T,2} # regarded as 2D
    panelmesh::PM
    submeshes::Vector{SM}
    size::NTuple{2,Int}
end

function CompositeMesh(panelmesh::PM, N) where {PM}
    T = eltype(panelmesh)
    submeshes = []
    for (i, p) in enumerate(panelmesh)
        intervals = AbstractMeshes.interval(panelmesh, i)
        DIM = length(intervals)
        origin = [intervals[j][1] for j in 1:DIM]
        lattice = diagm(DIM, DIM, [intervals[j][2] - intervals[j][1] for j in 1:DIM])
        #println(origin, lattice)
        cm = ChebMesh(origin, lattice, DIM, N)
        push!(submeshes, cm)
    end

    SM = typeof(submeshes[1])
    size = (length(submeshes[1]), length(panelmesh))
    return CompositeMesh{T,PM,SM}(panelmesh, submeshes, size)
end

Base.getindex(mesh::CompositeMesh, j, i) = mesh.submeshes[i][j]
Base.getindex(mesh::CompositeMesh, I::Int) = mesh[_ind2inds(size(mesh), I)...]

function AbstractMeshes.locate(mesh::CompositeMesh, x)
    pi = locate(mesh.panelmesh, x)
    si = locate(mesh.submeshes[pi], x)
    return _inds2ind(size(mesh), (si, pi))
end

AbstractMeshes.volume(mesh::CompositeMesh, j, i) = volume(mesh.submeshes[i], j)
AbstractMeshes.volume(mesh::CompositeMesh, I::Int) = volume(mesh, _ind2inds(size(mesh), I)...)
AbstractMeshes.volume(mesh::CompositeMesh) = volume(mesh.panelmesh)

AbstractMeshes.interval(mesh::CompositeMesh, j, i) = interval(mesh.submeshes[i], j)
AbstractMeshes.interval(mesh::CompositeMesh, I::Int) = interval(mesh, _ind2inds(size(mesh), I)...)

function AbstractMeshes.interp(data, mesh::CompositeMesh, x)
    pi = locate(mesh.panelmesh, x)
    dataview = reshape(view(data, :), size(mesh))
    slice = reshape(view(dataview, :, pi), size(mesh.submeshes[pi]))
    return interp(slice, mesh.submeshes[pi], x)
end

function AbstractMeshes.integrate(data, mesh::CompositeMesh)
    result = 0.0
    for (pi, p) in enumerate(mesh.panelmesh)
        dataview = reshape(view(data, :), size(mesh))
        slice = reshape(view(dataview, :, pi), size(mesh.submeshes[pi]))
        result += integrate(slice, mesh.submeshes[pi])
    end
    return result
end

end