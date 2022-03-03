module BaseMesh

using ..StaticArrays

export UniformMesh

abstract type AbstractMesh end

struct UniformMesh{DIM, N} <: AbstractMesh
    origin::SVector{DIM, Float64}
    latvec::SMatrix{DIM, DIM, Float64}
end

Base.length(mesh::UniformMesh{DIM, N}) where {DIM, N} = N
Base.size(mesh::UniformMesh{DIM, N}) where {DIM, N} = N^DIM
function Base.show(io::IO, mesh::UniformMesh)
    println("UniformMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end

# Base.view(mesh::UniformMesh, i::Int) = Base.view(mesh.mesh, :, i)

# # set is not allowed for meshs
function Base.getindex(mesh::UniformMesh{DIM, N}, inds...) where {DIM, N}
    pos = Vector(mesh.origin)
    for (ni, n) in enumerate(inds)
        pos = pos .+ mesh.latvec[ni, :] .* n ./ (length(mesh) + 1)
    end
    return pos
end
function Base.getindex(mesh::UniformMesh{DIM, N}, i::Int) where {DIM, N}
    inds = digits(i-1, base = N)
    pos = Vector(mesh.origin)
    for (ni, n) in enumerate(inds)
        pos = pos .+ mesh.latvec[ni, :] .* n ./ (N + 1)
    end
    return pos
end
Base.firstindex(mesh::UniformMesh) = 1
Base.lastindex(mesh::UniformMesh) = size(mesh)
# # iterator
Base.iterate(mesh::UniformMesh) = (mesh[1],1)
Base.iterate(mesh::UniformMesh, state) = (state>=size(mesh)) ? nothing : (mesh[state+1],state+1)

end
