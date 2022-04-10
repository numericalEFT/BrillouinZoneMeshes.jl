module BaseMesh

using ..StaticArrays

export UniformMesh

abstract type AbstractMesh end

struct UniformMesh{DIM, N} <: AbstractMesh
    origin::SVector{DIM, Float64}
    latvec::SMatrix{DIM, DIM, Float64}
    invlatvec::SMatrix{DIM, DIM, Float64}

    function UniformMesh{DIM, N}(origin, latvec) where {DIM, N}
        return new{DIM, N}(origin, latvec, inv(latvec))
    end
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
        pos = pos .+ mesh.latvec[:, ni] .* (n - 1 + 0.5) ./ (length(mesh))
    end
    return pos
end
function Base.getindex(mesh::UniformMesh{DIM, N}, i::Int) where {DIM, N}
    inds = digits(i-1, base = N, pad = DIM)
    pos = Vector(mesh.origin)
    for (ni, n) in enumerate(inds)
        pos = pos .+ mesh.latvec[:, ni] .* (n + 0.5) ./ (N)
    end
    return pos
end
Base.firstindex(mesh::UniformMesh) = 1
Base.lastindex(mesh::UniformMesh) = size(mesh)
# # iterator
Base.iterate(mesh::UniformMesh) = (mesh[1],1)
Base.iterate(mesh::UniformMesh, state) = (state>=size(mesh)) ? nothing : (mesh[state+1],state+1)

function Base.floor(mesh::UniformMesh{DIM, N}, x) where {DIM, N}
    # find index of nearest grid point to the point
    displacement = SVector{DIM, Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 1
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        if inds[i] < 1
            indexi = 1
        elseif inds[i] >= N
            indexi = N
        else
            indexi = floor(Int, inds[i])
        end
        # println("$(i):$(indexi)")
        indexall += (indexi - 1) * N ^ (i-1)
    end

    return indexall
end

end
