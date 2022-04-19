module BaseMesh

using ..StaticArrays
using ..LinearAlgebra

export UniformMesh, interp, integrate

abstract type AbstractMesh end

struct UniformMesh{DIM, N} <: AbstractMesh
    origin::SVector{DIM, Float64}
    latvec::SMatrix{DIM, DIM, Float64}
    invlatvec::SMatrix{DIM, DIM, Float64}

    function UniformMesh{DIM, N}(origin, latvec) where {DIM, N}
        return new{DIM, N}(origin, latvec, inv(latvec))
    end
end

Base.length(mesh::UniformMesh{DIM, N}) where {DIM, N} = N^DIM
Base.size(mesh::UniformMesh{DIM, N}) where {DIM, N} = NTuple{DIM, Int}(ones(Int, DIM) .* N)

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
Base.lastindex(mesh::UniformMesh) = length(mesh)
# # iterator
Base.iterate(mesh::UniformMesh) = (mesh[1],1)
Base.iterate(mesh::UniformMesh, state) = (state>=length(mesh)) ? nothing : (mesh[state+1],state+1)

_ind2inds(i::Int, N::Int, DIM::Int) = digits(i-1, base = N, pad = DIM) .+ 1
function _inds2ind(inds, N::Int)
    indexall = 1
    for i in 1:length(inds)
        indexall += (inds[i] - 1) * N ^ (i - 1)
    end
    return indexall
end
function _indfloor(x, N)
    if x < 1
        return 1
    elseif x >= N
        return N-1
    else
        return floor(Int, x)
    end
end

function Base.floor(mesh::UniformMesh{DIM, N}, x) where {DIM, N}
    # find index of nearest grid point to the point
    displacement = SVector{DIM, Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 0.5 .+ 2*eps(1.0*N)
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        indexall += (_indfloor(inds[i], N) - 1) * N ^ (i-1)
    end

    return indexall
end

function interp(data, mesh::UniformMesh, x)
    error("Not implemented!")
end

function interp(data::Matrix, mesh::UniformMesh{DIM, N}, x) where {DIM, N}
    @assert DIM == 2 "DIM should match dimension of data!"

    # find floor index and normalized x y
    ## find index of nearest grid point to the point
    displacement = SVector{DIM, Float64}(x) - mesh.origin
    xy = (mesh.invlatvec * displacement) .* N .+ 0.5 .+ 2*eps(N*1.0)
    xi, yi = _indfloor(xy[1], N), _indfloor(xy[2], N)

    return linear2D(data, xi, yi, xy[1], xy[2])
end

@inline function linear2D(data::Matrix, xi, yi, x, y)
    # accept data, floored index, normalized x and y, return linear interp
    # (xi, yi) should be [(1, 1) - size(data)], x and y normalized to the same scale as xi and yi
    dx0, dx1 = x - xi, xi - x + 1
    dy0, dy1 = y - yi, yi - y + 1

    d00, d01 = data[xi, yi], data[xi, yi+1]
    d10, d11 = data[xi+1, yi], data[xi+1, yi+1]

    g0 = d00 * dx1 + d10 * dx0
    g1 = d01 * dx1 + d11 * dx0

    gx = (g0 * dy1 + g1 * dy0) / (dx0 + dx1) / (dy0 + dy1)
    return gx
end

function integrate(data, mesh::UniformMesh{DIM, N}) where {DIM, N}
    area = abs(det(mesh.latvec))
    avg = sum(data) / length(data)
    return avg * area
end

end
