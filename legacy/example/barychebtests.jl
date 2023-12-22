using Test

module ChebMeshes

using BrillouinZoneMeshes
using BrillouinZoneMeshes.CompositeGrids
using BrillouinZoneMeshes.CompositeGrids.BaryChebTools

struct ChebMesh{T,DIM,S::Tuple} <: AbstractMesh{T,DIM}
    lattice::Matrix{T}
    inv_lattice::Matrix{T}
    cell_volume::T

    origin::SVector{DIM,T}
    size::NTuple{DIM,Int} # this should be (N,N,...)
    barychebs::Tuple

end

end

@testset "BaryCheb" begin

end