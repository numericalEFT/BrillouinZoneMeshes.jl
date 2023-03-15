
# this file is a part of BaseMesh.jl

#################### Abstract Orthogonal Mesh ###############

abstract type AbstractProdMesh{T,DIM} <: AbstractMesh{T,DIM} end
export AbstractProdMesh, DirectProdMesh, ProdMesh

_getgrid(mesh::AbstractProdMesh, I) = error("not implemented!")

function Base.getindex(mesh::AbstractProdMesh{T,DIM}, inds...) where {T,DIM}
    # i1, I = inds[1], AbstractMeshes._inds2ind(mesh.size, 1, inds[2:end]...)
    # seems generated function doesn't work here 
    # as compiler doesn't know mesh.size
    i1, I = inds[1], Base._sub2ind(size(mesh.mesh), inds[2:end]...)
    return SVector{DIM,T}([_getgrid(mesh, I)[i1], mesh.mesh[I]...])
end

function Base.getindex(mesh::AbstractProdMesh, I::Int)
    return Base.getindex(mesh, AbstractMeshes._ind2inds(mesh.size, I)...)
end

function AbstractMeshes.locate(mesh::AbstractProdMesh, x)
    I = AbstractMeshes.locate(mesh.mesh, x[2:end])
    i1 = AbstractMeshes.locate(_getgrid(mesh, I), x[1])
    return i1 + (I - 1) * size(mesh)[1]
end

function AbstractMeshes.volume(mesh::AbstractProdMesh)
    return sum(AbstractMeshes.volume(mesh.mesh, I) * AbstractMeshes.volume(_getgrid(mesh, I)) for I in 1:length(mesh.mesh))
end
function AbstractMeshes.volume(mesh::AbstractProdMesh, I::Int)
    i1, I2 = (I - 1) % size(mesh)[1] + 1, (I - 1) ÷ size(mesh)[1] + 1
    return AbstractMeshes.volume(mesh.mesh, I2) * AbstractMeshes.volume(_getgrid(mesh, I2), i1)
end

function AbstractMeshes.interval(mesh::AbstractProdMesh{T,DIM}, I::Int) where {T,DIM}
    inds = AbstractMeshes._ind2inds(size(mesh), I)
    J = Base._sub2ind(size(mesh)[2:end], inds[2:end]...)
    x1, x2 = AbstractMeshes.interval(_getgrid(mesh, J), inds[1])
    result = MMatrix{DIM,2,T,DIM * 2}(zeros(DIM, 2))
    if DIM == 2
        result[1, :] .= (x1, x2)
        result[2, :] .= AbstractMeshes.interval(mesh.mesh, J)
    else
        result[1, :] .= (x1, x2)
        result[2:end, :] .= AbstractMeshes.interval(mesh.mesh, J)
    end
    return SMatrix{DIM,2,T,DIM * 2}(result)
end

# DirectProdMesh is mesh from direct product grid×mesh
struct DirectProdMesh{T,DIM,MT,GT<:AbstractGrid} <: AbstractProdMesh{T,DIM}
    mesh::MT
    grid::GT
    size::NTuple{DIM,Int}
end

function DirectProdMesh(grid::GT, mesh::AbstractMesh{T,DIM}) where {T,DIM,GT}
    MT = typeof(mesh)
    # N of element in grids should match length of mesh
    msize = (length(grid), size(mesh)...)
    return DirectProdMesh{T,DIM + 1,MT,GT}(mesh, grid, msize)
end

function DirectProdMesh(grid::GT, mesh::MT) where {MT<:AbstractGrid,GT}
    msize = (length(grid), length(mesh))
    return DirectProdMesh{eltype(MT),2,MT,GT}(mesh, grid, msize)
end

function DirectProdMesh(grids...)
    # this method allow constructor like DirectProdMesh(r, theta, phi)
    # where all arguments are AbstractGrid
    if length(grids) <= 2
        error("this method shouldn't be called for less than 2 grids!")
    else
        return DirectProdMesh(grids[1], DirectProdMesh(grids[2:end]...))
    end
end

_getgrid(mesh::DirectProdMesh, I) = mesh.grid

# ProdMesh is mesh from grids×mesh
# ProdMesh is generalized version that allows different grid for different meshpoints
# equal length first. 
# ProdMesh without equal length makes linear indexing difficult
struct ProdMesh{T,DIM,MT,GT<:AbstractGrid} <: AbstractProdMesh{T,DIM}
    # composite mesh constructed upon a mesh and a set of grids
    # the dimension represented by the grids becomes the first dimension
    # while other dimensions are represented by mesh
    # the mesh point at (i,j,k...) will be (mesh.grids[j,k...][i], mesh.mesh[j,k...]...)
    mesh::MT
    grids::Vector{GT}
    size::NTuple{DIM,Int}
end

function ProdMesh(grids::Vector{GT}, mesh::AbstractMesh{T,DIM}) where {T,DIM,GT}
    MT = typeof(mesh)
    # N of element in grids should match length of mesh
    @assert length(mesh) == length(grids)
    # length of all grids should be the same
    @assert length.(grids) == ones(length(grids)) .* length(grids[1])
    msize = (length(grids[1]), size(mesh)...)
    return ProdMesh{T,DIM + 1,MT,GT}(mesh, grids, msize)
end

function ProdMesh(grids::Vector{GT}, mesh::MT) where {MT<:AbstractGrid,GT}
    @assert length(mesh) == length(grids)
    @assert length.(grids) == ones(length(grids)) .* length(grids[1])
    msize = (length(grids[1]), length(mesh))
    return ProdMesh{eltype(MT),2,MT,GT}(mesh, grids, msize)
end

_getgrid(mesh::ProdMesh, I) = mesh.grids[I]
