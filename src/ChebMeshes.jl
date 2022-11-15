
# this file is included in BaseMesh.jl

struct ChebMesh{T,DIM,N} <: AbstractMesh{T,DIM}
    lattice::Matrix{T}
    inv_lattice::Matrix{T}
    cell_volume::T

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int}
    barycheb::BaryCheb1D{N}

    function ChebMesh(origin, lattice, DIM, bc::BaryCheb1D{N}) where {N}
        T = eltype(origin)
        inv_lattice = inv(lattice)
        cell_volume = abs(det(lattice))
        size = NTuple{DIM,Int}(N for i in 1:N)

        return new{T,DIM,N}(lattice, inv_lattice, cell_volume, origin, size, bc)
    end
end

# following constructors take generated barycheb to avoid waste
function ChebMesh(origin, lattice, DIM, N::Int)
    return ChebMesh(origin, lattice, DIM, BaryCheb1D(N))
end
function ChebMesh(origin, lattice, cm::ChebMesh{T,DIM,N}) where {T,DIM,N}
    return ChebMesh(origin, lattice, DIM, cm.barycheb)
end

# handle latvec with HasLatAndInv
AbstractMeshes.MeshDomain(::Type{<:ChebMesh}) = OnLattice()

# Base.length(mesh::ChebMesh) = prod(mesh.size)
# Base.size(mesh::ChebMesh) = mesh.size
# Base.size(mesh::ChebMesh, I) = mesh.size[I]

function Base.getindex(mesh::ChebMesh{T,DIM,N}, inds...) where {T,DIM,N}
    a = SVector{DIM,Float64}(mesh.barycheb[n] for n in inds)
    return mesh.origin + AbstractMeshes.frac_to_cart(mesh, (a .+ 1.0) ./ 2.0)
end
Base.getindex(mesh::ChebMesh, I::Int) = mesh[AbstractMeshes._ind2inds(mesh.size, I)...]

function AbstractMeshes.interp(data, mesh::ChebMesh{T,DIM,N}, x) where {T,DIM,N}
    # translate into dimensionless
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = AbstractMeshes.cart_to_frac(mesh, displacement) .* 2.0 .- 1.0
    return interpND(data, mesh.barycheb, xs)
end

function AbstractMeshes.integrate(data, mesh::ChebMesh{T,DIM,N}) where {T,DIM,N}
    return integrateND(data, mesh.barycheb, DIM) * mesh.cell_volume / 2^DIM
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

function AbstractMeshes.locate(mesh::ChebMesh{T,DIM,N}, x) where {T,DIM,N}
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = AbstractMeshes.cart_to_frac(displacement) .* 2.0 .- 1.0
    inds = NTuple{DIM,Int}(locate1d(mesh.barycheb, xs[i]) for i in 1:DIM)
    return AbstractMeshes._inds2ind(mesh.size, inds)
end

AbstractMeshes.volume(mesh::ChebMesh) = mesh.cell_volume
function AbstractMeshes.volume(mesh::ChebMesh{T,DIM,N}, i) where {T,DIM,N}
    inds = AbstractMeshes._ind2inds(mesh.size, i)
    return reduce(*, volume1d(mesh.barycheb, inds[j]) for j in 1:DIM) * volume(mesh) / 2^DIM
end

