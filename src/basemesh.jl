module BaseMesh

using ..StaticArrays
using ..LinearAlgebra

using ..BaryCheb

export UniformMesh, BaryChebMesh

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

function interp(data, mesh::UniformMesh{DIM, N}, x) where {DIM, N}
    error("Not implemented!")
end

function interp(data, mesh::UniformMesh{2, N}, x) where {N}
    # find floor index and normalized x y
    displacement = SVector{2, Float64}(x) - mesh.origin
    xy = (mesh.invlatvec * displacement) .* N .+ 0.5 .+ 2*eps(N*1.0)
    xi, yi = _indfloor(xy[1], N), _indfloor(xy[2], N)

    return linear2D(data, xi, yi, xy..., N)
end

@inline function linear2D(data, xi, yi, x, y, N)
    # accept data, floored index, normalized x and y, return linear interp
    # (xi, yi) should be [(1, 1) - size(data)], x and y normalized to the same scale as xi and yi
    xd, yd= x-xi, y-yi

    c0 = data[xi + (yi-1) * N] * (1-xd) + data[xi+1 + (yi-1) * N] * xd
    c1 = data[xi + (yi) * N] * (1-xd) + data[xi+1 + (yi) * N] * xd

    return c0 * (1-yd) + c1 * yd
end

function interp(data, mesh::UniformMesh{3, N}, x) where {T, N}
    # find floor index and normalized x y z
    displacement = SVector{3, Float64}(x) - mesh.origin
    xyz = (mesh.invlatvec * displacement) .* N .+ 0.5 .+ 2*eps(N*1.0)
    xi, yi, zi = _indfloor(xyz[1], N), _indfloor(xyz[2], N), _indfloor(xyz[3], N)

    return linear3D(data, xi, yi, zi, xyz..., N)
end

@inline function linear3D(data, xi, yi, zi, x, y, z, N) where {T}
    xd, yd, zd = x-xi, y-yi, z-zi

    c00 = data[xi + (yi-1)*N + (zi-1)*N^2] * (1-xd) + data[xi + (yi-1)*N + (zi-1)*N^2] * xd
    c01 = data[xi + (yi-1)*N + (zi)*N^2] * (1-xd) + data[xi + (yi-1)*N + (zi)*N^2] * xd
    c10 = data[xi + (yi)*N + (zi-1)*N^2] * (1-xd) + data[xi + (yi)*N + (zi-1)*N^2] * xd
    c11 = data[xi + (yi)*N + (zi)*N^2] * (1-xd) + data[xi + (yi)*N + (zi)*N^2] * xd

    c0 = c00 * (1-yd) + c10 * yd
    c1 = c01 * (1-yd) + c11 * yd

    return c0 * (1-zd) + c1 * zd
end

function integrate(data, mesh::UniformMesh{DIM, N}) where {DIM, N}
    area = abs(det(mesh.latvec))
    avg = sum(data) / length(data)
    return avg * area
end

struct BaryChebMesh{DIM, N} <: AbstractMesh
    origin::SVector{DIM, Float64}
    latvec::SMatrix{DIM, DIM, Float64}
    invlatvec::SMatrix{DIM, DIM, Float64}
    barycheb::BaryCheb1D{N} # grid in barycheb is [-1, 1]

    function BaryChebMesh{DIM, N}(origin, latvec, barycheb::BaryCheb1D{N}) where {DIM, N}
        return new{DIM, N}(origin, latvec, inv(latvec), barycheb)
    end
end

function BaryChebMesh(origin, latvec, DIM, N)
    barycheb = BaryCheb1D(N)
    return BaryChebMesh{DIM, N}(origin, latvec, barycheb)
end

Base.length(mesh::BaryChebMesh{DIM, N}) where {DIM, N} = N^DIM
Base.size(mesh::BaryChebMesh{DIM, N}) where {DIM, N} = NTuple{DIM, Int}(ones(Int, DIM) .* N)
function Base.show(io::IO, mesh::BaryChebMesh)
    println("BaryChebMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end
function Base.getindex(mesh::BaryChebMesh{DIM, N}, inds...) where {DIM, N}
    pos = Vector(mesh.origin)
    for (ni, n) in enumerate(inds)
        pos = pos .+ mesh.latvec[:, ni] .* (mesh.barycheb[n] + 1.0) ./ 2.0
    end
    return pos
end
function Base.getindex(mesh::BaryChebMesh{DIM, N}, i::Int) where {DIM, N}
    inds = digits(i-1, base = N, pad = DIM)
    pos = Vector(mesh.origin)
    for (ni, n) in enumerate(inds)
        pos = pos .+ mesh.latvec[:, ni] .* (mesh.barycheb[n+1] + 1.0) ./ 2.0
    end
    return pos
end
Base.firstindex(mesh::BaryChebMesh) = 1
Base.lastindex(mesh::BaryChebMesh) = length(mesh)
# # iterator
Base.iterate(mesh::BaryChebMesh) = (mesh[1],1)
Base.iterate(mesh::BaryChebMesh, state) = (state>=length(mesh)) ? nothing : (mesh[state+1],state+1)

function interp(data, mesh::BaryChebMesh{DIM, N}, x) where {DIM, N}
    # translate into dimensionless
    displacement = SVector{DIM, Float64}(x) - mesh.origin
    xs = (mesh.invlatvec * displacement) .* 2.0 .- 1.0
    return interpND(data, mesh.barycheb, xs)
end

function integrate(data, mesh::BaryChebMesh{DIM, N}) where {DIM, N}
    area = abs(det(mesh.latvec))
    return integrateND(data, mesh.barycheb, DIM) * area / 2^DIM
end

end
  
