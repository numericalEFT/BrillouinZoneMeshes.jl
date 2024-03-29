module PointSymmetry
using spglib_jll
using LinearAlgebra
using Printf
using StaticArrays
using Spglib
export get_ir_reciprocal_mesh
export standardize_cell


# Deduce the value type (floating-point type for storing plain an static data
# in Model and PlaneWaveBasis) from e.g. an interval or a dual type.
value_type(T) = T

# Frequently-used array types
const Mat3{T} = SMatrix{3,3,T,9} where {T}
const Vec3{T} = SVector{3,T} where {T}
const AbstractArray3{T} = AbstractArray{T,3}

function _makeVec3(kcoord::AbstractVector{T}) where {T}
    if length(kcoord) == 3
        return Vec3{T}(kcoord)
    elseif length(kcoord) == 2
        return Vec3{T}(kcoord[1], kcoord[2], 0)
    elseif length(kcoord) == 1
        return Vec3{T}(kcoord[1], 0, 0)
    else
        error("illegal kcoord dimension!")
    end
end

function _make3D(lattice::AbstractMatrix{T}, position::AbstractVector) where {T}
    @assert size(lattice, 1) == size(lattice, 2)
    DIM = size(lattice, 1)
    _lattice = zeros(T, 3, 3)
    for i in 1:DIM
        _lattice[i, 1:DIM] = lattice[i, 1:DIM]
    end
    for i in DIM+1:3
        _lattice[i, i] = 1
    end
    _position = [zeros(T, 3) for i in 1:length(position)]
    for ai in eachindex(_position)
        _position[ai][1:DIM] = position[ai][1:DIM]
    end
    return _lattice, _position
end

include("SymOp.jl")
include("spglib.jl")
include("symmetry.jl")

end