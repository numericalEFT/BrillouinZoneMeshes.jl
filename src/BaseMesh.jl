module BaseMesh

using ..StaticArrays
using ..LinearAlgebra
using ..CompositeGrids

using ..BaryCheb
using ..AbstractMeshes
using ..Model
using ..Model: get_latvec


export UniformMesh, BaryChebMesh, CenteredMesh, EdgedMesh, UMesh, CompositeMesh

############## Abstract Uniform Mesh #################
abstract type AbstractUniformMesh{T,DIM} <: AbstractMesh{T,DIM} end
export AbstractUniformMesh
export inv_lattice_vector, lattice_vector, cell_volume

Base.length(mesh::AbstractUniformMesh) = prod(mesh.size)
Base.size(mesh::AbstractUniformMesh) = mesh.size
Base.size(mesh::AbstractUniformMesh, I) = mesh.size[I]

function AbstractMeshes.fractional_coordinates(mesh::AbstractUniformMesh{T,DIM}, I::Int) where {T,DIM}
    n = SVector{DIM,Int}(AbstractMeshes._ind2inds(mesh.size, I))
    return (n .- 1 .+ mesh.shift) ./ mesh.size
end

function AbstractMeshes.fractional_coordinates(mesh::AbstractUniformMesh{T,DIM}, x::AbstractVector) where {T,DIM}
    displacement = SVector{DIM,T}(x)
    return (inv_lattice_vector(mesh) * displacement)
end

function Base.getindex(mesh::AbstractUniformMesh{T,DIM}, inds...) where {T,DIM}
    n = SVector{DIM,Int}(inds)
    return mesh.origin + lattice_vector(mesh) * ((n .- 1 .+ mesh.shift) ./ mesh.size)
end

function Base.getindex(mesh::AbstractUniformMesh, I::Int)
    return Base.getindex(mesh, AbstractMeshes._ind2inds(mesh.size, I)...)
end

"""
    function AbstractMeshes.locate(mesh::AbstractUniformMesh{T,DIM}, x) where {T,DIM}

locate mesh point in mesh that is nearest to x. Useful for Monte-Carlo algorithm.
Could also be used for zeroth order interpolation.

# Parameters
- `mesh`: aimed mesh
- `x`: cartesian pos to locate
"""
function AbstractMeshes.locate(mesh::AbstractUniformMesh{T,DIM}, x) where {T,DIM}
    # find index of nearest grid point to the point
    svx = SVector{DIM,T}(x)
    inds = fractional_coordinates(mesh, svx - mesh.origin) .* mesh.size .+ 1.5 .- mesh.shift .+ 2 .* eps.(T.(mesh.size))
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    factor = 1
    # indexall += (_indfloor(inds[1], mesh.size[1]; edgeshift=0) - 1) * factor
    indexall += (cycling_floor(inds[1], mesh.size[1]) - 1) * factor
    for i in 2:DIM
        factor *= mesh.size[i-1]
        indexall += (cycling_floor(inds[i], mesh.size[i]) - 1) * factor
    end

    return indexall
end

"""
    function AbstractMeshes.volume(mesh::AbstractUniformMesh, i)

volume represented by mesh point i. When i is omitted return volume of the whole mesh. 
For M-P mesh it's always volume(mesh)/length(mesh), but for others things are more complecated.
Here we assume periodic boundary condition so for all case it's the same.

# Parameters:
- `mesh`: mesh
- `i`: index of mesh point, if ommited return volume of whole mesh
"""
AbstractMeshes.volume(mesh::AbstractUniformMesh) = cell_volume(mesh)
AbstractMeshes.volume(mesh::AbstractUniformMesh, i) = cell_volume(mesh) / length(mesh)

############## general purposed uniform mesh #################
struct UMesh{T,DIM} <: AbstractUniformMesh{T,DIM}
    lattice::Matrix{T}
    inv_lattice::Matrix{T}
    cell_volume::T

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    shift::SVector{DIM,Rational}
end

# UMesh(;
#     br::Brillouin{T,DIM},
#     origin::Real,
#     size,
#     shift::Real) where {T,DIM} = UMesh{T,DIM}(
#     br.recip_lattice,
#     br.inv_recip_lattice,
#     br.recip_cell_volume,
#     SVector{DIM,T}(br.recip_lattice * ones(T, DIM) .* origin),
#     size,
#     SVector{DIM,Rational}(shift .* ones(Int, DIM))
# )

UMesh(;
    br::Brillouin{T,DIM},
    origin,
    size,
    shift) where {T,DIM} = UMesh{T,DIM}(
    br.recip_lattice,
    br.inv_recip_lattice,
    br.recip_cell_volume,
    SVector{DIM,T}(br.recip_lattice * origin),
    size,
    SVector{DIM,Rational}(shift)
)

lattice_vector(mesh::UMesh) = mesh.lattice
inv_lattice_vector(mesh::UMesh) = mesh.inv_lattice
lattice_vector(mesh::UMesh, i::Int) = get_latvec(mesh.lattice, i)
inv_lattice_vector(mesh::UMesh, i::Int) = get_latvec(mesh.inv_lattice, i)
cell_volume(mesh::UMesh) = mesh.cell_volume

function Base.show(io::IO, mesh::UMesh)
    println("UMesh with $(length(mesh)) mesh points")
end

function cycling_floor(I, N)
    ifloor = (floor(Int, I) + N) % N
    if ifloor == 0
        return N
    else
        return ifloor
    end
end


# equal length first. 
# CompositeMesh without equal length makes linear indexing difficult
struct CompositeMesh{T,DIM,MT,GT<:AbstractGrid} <: AbstractMesh{T,DIM}
    # composite mesh constructed upon a mesh and a set of grids
    # the dimension represented by the grids becomes the first dimension
    # while other dimensions are represented by mesh
    # the mesh point at (i,j,k...) will be (mesh.grids[j,k...][i], mesh.mesh[j,k...]...)
    mesh::MT
    grids::Vector{GT}
    size::NTuple{DIM,Int}
end

function CompositeMesh(mesh::AbstractMesh{T,DIM}, grids::Vector{GT}) where {T,DIM,GT}
    MT = typeof(mesh)
    # N of element in grids should match length of mesh
    @assert length(mesh) == length(grids)
    # length of all grids should be the same
    @assert length.(grids) == ones(length(grids)) .* length(grids[1])
    msize = (length(grids[1]), size(mesh)...)
    return CompositeMesh{T,DIM + 1,MT,GT}(mesh, grids, msize)
end

function CompositeMesh(mesh::MT, grids::Vector{GT}) where {MT<:AbstractGrid,GT}
    @assert length(mesh) == length(grids)
    @assert length.(grids) == ones(length(grids)) .* length(grids[1])
    msize = (length(grids[1]), length(mesh))
    return CompositeMesh{eltype(MT),2,MT,GT}(mesh, grids, msize)
end

Base.length(mesh::CompositeMesh) = prod(mesh.size)
Base.size(mesh::CompositeMesh) = mesh.size
Base.size(mesh::CompositeMesh, I::Int) = mesh.size[I]

function Base.getindex(mesh::CompositeMesh{T,DIM,MT,GT}, inds...) where {T,DIM,MT,GT}
    # i1, I = inds[1], AbstractMeshes._inds2ind(mesh.size, 1, inds[2:end]...)
    # seems generated function doesn't work here 
    # as compiler doesn't know mesh.size
    i1, I = inds[1], Base._sub2ind(size(mesh.mesh), inds[2:end]...)
    return SVector{DIM,T}([mesh.grids[I][i1], mesh.mesh[I]...])
end

function Base.getindex(mesh::CompositeMesh, I::Int)
    return Base.getindex(mesh, AbstractMeshes._ind2inds(mesh.size, I)...)
end

function AbstractMeshes.locate(mesh::CompositeMesh, x)
    I = AbstractMeshes.locate(mesh.mesh, x[2:end])
    i1 = AbstractMeshes.locate(mesh.grids[I], x[1])
    return i1 + (I - 1) * size(mesh)[1]
end

function AbstractMeshes.volume(mesh::CompositeMesh)
    return sum(AbstractMeshes.volume(mesh.mesh, I) * AbstractMeshes.volume(mesh.grids[I]) for I in 1:length(mesh.mesh))
end
function AbstractMeshes.volume(mesh::CompositeMesh, I::Int)
    i1, I2 = (I - 1) % size(mesh)[1] + 1, (I - 1) ÷ size(mesh)[1] + 1
    return AbstractMeshes.volume(mesh.mesh, I2) * AbstractMeshes.volume(mesh.grids[I2], i1)
end

#####################################
# LEGACY CODE BELOW
#####################################

abstract type EqualLengthMesh{DIM,N} <: AbstractMesh{Float64,DIM} end

abstract type MeshType end
struct CenteredMesh <: MeshType end # Monkhorst-Pack mesh, take center points instead of left-bottom
struct EdgedMesh <: MeshType end # Γ=(0,0) centered, take (0,0) as mesh point

# TODO: support (N1, N2, N3)
struct UniformMesh{DIM,N,MT,DIMSQ} <: EqualLengthMesh{DIM,N}
    # DIMSQ should be provided explicitly in SMatrix parameters
    # otherwise it cause type instability
    origin::SVector{DIM,Float64}
    latvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    invlatvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    # dims could be here as field element

    function UniformMesh{DIM,N,MT}(origin, latvec) where {DIM,N,MT<:MeshType}
        return new{DIM,N,MT,DIM^2}(origin, latvec, inv(latvec))
    end
end

function UniformMesh{DIM,N}(origin, latvec) where {DIM,N}
    # Centered Mesh is the default
    return UniformMesh{DIM,N,CenteredMesh}(origin, latvec)
end

Base.length(mesh::UniformMesh{DIM,N}) where {DIM,N} = N^DIM
Base.size(mesh::UniformMesh{DIM,N}) where {DIM,N} = NTuple{DIM,Int}(ones(Int, DIM) .* N)

function Base.show(io::IO, mesh::UniformMesh)
    println("UniformMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end

# Base.view(mesh::UniformMesh, i::Int) = Base.view(mesh.mesh, :, i)

# # set is not allowed for meshs
meshshift(::Type) = error("not implement!")
meshshift(::Type{<:CenteredMesh}) = 0.5
meshshift(::Type{<:EdgedMesh}) = 0.0

function Base.getindex(mesh::UniformMesh{DIM,N,MT}, inds...) where {DIM,N,MT}
    # pos = Vector(mesh.origin)
    # for (ni, n) in enumerate(inds)
    #     pos = pos .+ mesh.latvec[:, ni] .* (n - 1 + meshshift(MT)) ./ (N)
    # end
    # return pos
    n = SVector{DIM,Int}(inds)
    return mesh.origin + mesh.latvec * ((n .- 1 .+ meshshift(MT)) ./ N)
end

function Base.getindex(mesh::UniformMesh{DIM,N,MT}, i::Int) where {DIM,N,MT}
    return Base.getindex(mesh, _ind2inds(mesh, i)...)
end
Base.firstindex(mesh::UniformMesh) = 1
Base.lastindex(mesh::UniformMesh) = length(mesh)
# # iterator
Base.iterate(mesh::UniformMesh) = (mesh[1], 1)
Base.iterate(mesh::UniformMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

# _ind2inds(i::Int, N::Int, DIM::Int) = digits(i - 1, base=N, pad=DIM) .+ 1

# function _inds2ind(inds, N::Int)
#     indexall = 1
#     for i in 1:length(inds)
#         indexall += (inds[i] - 1) * N^(i - 1)
#     end
#     return indexall
# end

@generated function _inds2ind(umesh::EqualLengthMesh{DIM,N}, I) where {DIM,N}
    ex = :(I[DIM] - 1)
    for i = (DIM-1):-1:1
        ex = :(I[$i] - 1 + N * $ex)
    end
    return :($ex + 1)
end

@generated function _ind2inds(umesh::EqualLengthMesh{DIM,N}, I::Int) where {DIM,N}
    inds, quotient = :((I - 1) % N + 1), :((I - 1) ÷ N)
    for i = 2:DIM-1
        inds, quotient = :($inds..., $quotient % N + 1), :($quotient ÷ N)
    end
    inds = :($inds..., $quotient + 1)
    return :(SVector{DIM,Int}($inds))
end

function _indfloor(x, N; edgeshift=1)

    # edgeshift = 1 by default in floor function so that end point return N-1
    # edgeshift = 0 in locate function
    if x < 1
        return 1
    elseif x >= N
        return N - edgeshift
    else
        return floor(Int, x)
    end
end

function Base.floor(mesh::UniformMesh{DIM,N}, x) where {DIM,N}
    # find index of nearest grid point to the point
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 0.5 .+ 2 * eps(1.0 * N)
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        indexall += (_indfloor(inds[i], N) - 1) * N^(i - 1)
    end

    return indexall
end

function AbstractMeshes.locate(mesh::UniformMesh{DIM,N,MT}, x) where {DIM,N,MT}
    # find index of nearest grid point to the point
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    # println(displacement)
    inds = (mesh.invlatvec * displacement) .* (N) .+ 1.5 .- meshshift(MT) .+ 2 * eps(1.0 * N)
    indexall = 1
    # println((mesh.invlatvec * displacement))
    # println(inds)
    for i in 1:DIM
        indexall += (_indfloor(inds[i], N; edgeshift=0) - 1) * N^(i - 1)
    end
    return indexall
end

AbstractMeshes.volume(mesh::UniformMesh) = abs(det(mesh.latvec))
function AbstractMeshes.volume(mesh::UniformMesh{DIM,N,MT}, i) where {DIM,N,MT}
    volume(MT, mesh, i)
end

function AbstractMeshes.volume(::Type{<:CenteredMesh}, mesh, i)
    # for uniform centered mesh, all mesh points share the same volume
    return abs(det(mesh.latvec)) / length(mesh)
end

function AbstractMeshes.volume(::Type{<:EdgedMesh}, mesh::UniformMesh{DIM,N,MT}, i) where {DIM,N,MT}
    inds = _ind2inds(mesh, i)
    n1, nend = count(i -> i == 1, inds), count(i -> i == N, inds)
    cellarea = 2^(DIM - n1 - nend) * 3^nend / 2^DIM
    return cellarea / length(mesh) * volume(mesh)
end

function interp(data, mesh::UniformMesh{DIM,N}, x) where {DIM,N}
    error("Not implemented!")
end

function interp(data, mesh::UniformMesh{2,N,MT}, x) where {N,MT}
    # find floor index and normalized x y
    displacement = SVector{2,Float64}(x) - mesh.origin
    xy = (mesh.invlatvec * displacement) .* N .+ (1 - meshshift(MT)) .* (1 + 2 * eps(N * 1.0))
    xi, yi = _indfloor(xy[1], N), _indfloor(xy[2], N)

    return linear2D(data, xi, yi, xy..., N)
end

@inline function linear2D(data, xi, yi, x, y, N)
    # accept data, floored index, normalized x and y, return linear interp
    # (xi, yi) should be [(1, 1) - size(data)], x and y normalized to the same scale as xi and yi
    xd, yd = x - xi, y - yi

    c0 = data[xi+(yi-1)*N] * (1 - xd) + data[xi+1+(yi-1)*N] * xd
    c1 = data[xi+(yi)*N] * (1 - xd) + data[xi+1+(yi)*N] * xd

    return c0 * (1 - yd) + c1 * yd
end

function interp(data, mesh::UniformMesh{3,N,MT}, x) where {T,N,MT}
    # find floor index and normalized x y z
    displacement = SVector{3,Float64}(x) - mesh.origin
    xyz = (mesh.invlatvec * displacement) .* N .+ (1 - meshshift(MT)) .* (1 + 4 * eps(N * 1.0))
    xi, yi, zi = _indfloor(xyz[1], N), _indfloor(xyz[2], N), _indfloor(xyz[3], N)

    return linear3D(data, xi, yi, zi, xyz..., N)
end

@inline function linear3D(data, xi, yi, zi, x, y, z, N) where {T}
    xd, yd, zd = x - xi, y - yi, z - zi

    c00 = data[xi+(yi-1)*N+(zi-1)*N^2] * (1 - xd) + data[xi+(yi-1)*N+(zi-1)*N^2] * xd
    c01 = data[xi+(yi-1)*N+(zi)*N^2] * (1 - xd) + data[xi+(yi-1)*N+(zi)*N^2] * xd
    c10 = data[xi+(yi)*N+(zi-1)*N^2] * (1 - xd) + data[xi+(yi)*N+(zi-1)*N^2] * xd
    c11 = data[xi+(yi)*N+(zi)*N^2] * (1 - xd) + data[xi+(yi)*N+(zi)*N^2] * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    return c0 * (1 - zd) + c1 * zd
end

function integrate(data, mesh::UniformMesh{DIM,N,CenteredMesh}) where {DIM,N}
    area = abs(det(mesh.latvec))
    avg = sum(data) / length(data)
    return avg * area
end

function integrate(data, mesh::UniformMesh{DIM,N,EdgedMesh}) where {DIM,N}
    area = abs(det(mesh.latvec))
    avg = 0.0
    for i in 1:length(data)
        inds = _ind2inds(mesh, i)
        n1, nend = count(i -> i == 1, inds), count(i -> i == N, inds)
        avg += 2^(DIM - n1 - nend) * 3^nend / 2^DIM * data[i]
    end
    avg = avg / length(data)
    return avg * area
end

struct BaryChebMesh{DIM,N,DIMSQ} <: EqualLengthMesh{DIM,N}
    origin::SVector{DIM,Float64}
    latvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    invlatvec::SMatrix{DIM,DIM,Float64,DIMSQ}
    barycheb::BaryCheb1D{N} # grid in barycheb is [-1, 1]

    function BaryChebMesh{DIM,N}(origin, latvec, barycheb::BaryCheb1D{N}) where {DIM,N}
        return new{DIM,N,DIM^2}(origin, latvec, inv(latvec), barycheb)
    end
end

function BaryChebMesh(origin, latvec, DIM, N)
    barycheb = BaryCheb1D(N)
    return BaryChebMesh{DIM,N}(origin, latvec, barycheb)
end

function BaryChebMesh(origin, latvec, bcmesh::BaryChebMesh{DIM,N}) where {DIM,N}
    # this constructor borrows barycheb from generated bcmesh
    barycheb = bcmesh.barycheb
    return BaryChebMesh{DIM,N}(origin, latvec, barycheb)
end

Base.length(mesh::BaryChebMesh{DIM,N}) where {DIM,N} = N^DIM
Base.size(mesh::BaryChebMesh{DIM,N}) where {DIM,N} = NTuple{DIM,Int}(ones(Int, DIM) .* N)
function Base.show(io::IO, mesh::BaryChebMesh)
    println("BaryChebMesh Grid:")
    for (pi, p) in enumerate(mesh)
        println(p)
    end
end
function Base.getindex(mesh::BaryChebMesh{DIM,N}, inds...) where {DIM,N}
    # pos = Vector(mesh.origin)
    # for (ni, n) in enumerate(inds)
    #     pos = pos .+ mesh.latvec[:, ni] .* (mesh.barycheb[n] + 1.0) ./ 2.0
    # end
    # return pos
    a = SVector{DIM,Float64}(mesh.barycheb[n] for n in inds)
    return mesh.origin + mesh.latvec * (a .+ 1.0) ./ 2.0
end
function Base.getindex(mesh::BaryChebMesh{DIM,N}, i::Int) where {DIM,N}
    return Base.getindex(mesh, _ind2inds(mesh, i)...)
end

Base.firstindex(mesh::BaryChebMesh) = 1
Base.lastindex(mesh::BaryChebMesh) = length(mesh)
# # iterator
Base.iterate(mesh::BaryChebMesh) = (mesh[1], 1)
Base.iterate(mesh::BaryChebMesh, state) = (state >= length(mesh)) ? nothing : (mesh[state+1], state + 1)

function interp(data, mesh::BaryChebMesh{DIM,N}, x) where {DIM,N}
    # translate into dimensionless
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = (mesh.invlatvec * displacement) .* 2.0 .- 1.0
    return interpND(data, mesh.barycheb, xs)
end

function integrate(data, mesh::BaryChebMesh{DIM,N}) where {DIM,N}
    area = abs(det(mesh.latvec))
    return integrateND(data, mesh.barycheb, DIM) * area / 2^DIM
end

function locate1d(g::BaryCheb1D, x)
    grid = g.x
    if x <= grid[1]
        return 1
    elseif x >= grid[end]
        if length(grid) != 1
            return length(grid)
        end
    end

    i2 = searchsortedfirst(grid, x)
    i1 = i2 - 1
    return abs(grid[i1] - x) < abs(grid[i2] - x) ? i1 : i2
end

function volume1d(g::BaryCheb1D, i)
    grid = g.x
    if i != 1 && i != length(grid)
        return (grid[i+1] - grid[i-1]) / 2
    elseif i == 1
        return (grid[i+1] + grid[i]) / 2 - (-1)
    else
        return 1 - (grid[i] + grid[i-1]) / 2
    end
end

function AbstractMeshes.locate(mesh::BaryChebMesh{DIM,N}, x) where {DIM,N}
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = (mesh.invlatvec * displacement) .* 2.0 .- 1.0
    inds = [locate1d(mesh.barycheb, xs[i]) for i in 1:DIM]
    return _inds2ind(mesh, inds)
end

AbstractMeshes.volume(mesh::BaryChebMesh) = abs(det(mesh.latvec))
function AbstractMeshes.volume(mesh::BaryChebMesh{DIM,N}, i) where {DIM,N}
    inds = _ind2inds(mesh, i)
    return reduce(*, volume1d(mesh.barycheb, inds[j]) for j in 1:DIM) * volume(mesh) / 2^DIM
end

end

