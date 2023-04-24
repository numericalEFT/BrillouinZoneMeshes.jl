module BrillouinZoneMeshes

using AbstractTrees
using StaticArrays
using Statistics
using LinearAlgebra
using CompositeGrids
using Printf
using Roots
# using CoordinateTransformations

# Write your package code here.
BaryCheb = CompositeGrids.BaryChebTools
# export BaryCheb

include("AbstractMeshes.jl")
using .AbstractMeshes
export AbstractMeshes, AbstractMesh
export fractional_coordinates

include("printing.jl")

include("symmetry/PointSymmetry.jl")
export PointSymmetry

include("Cells.jl")
using .Cells
export Cells

include("BaseMesh.jl")
using .BaseMesh
export BaseMesh
export UniformMesh, BaryChebMesh, CenteredMesh, EdgedMesh, AbstractMesh# , locate, volume
export inv_lattice_vector, lattice_vector, cell_volume
export AbstractUniformMesh

include("CompositeMeshes.jl")
using .CompositeMeshes
# export CompositeMeshes

# include("TreeMeshes.jl")
# using .TreeMeshes
# export TreeMeshes
# export GridNode, TreeGrid, uniformtreegrid, treegridfromdensity, efficiency# , interp, integrate

include("MeshMaps.jl")
using .MeshMaps
export MeshMaps
export SymMap, MappedData, MeshMap, ReducedBZMesh

include("BZMeshes.jl")
using .BZMeshes
export BZMeshes
export UniformBZMesh


using Requires



function __init__()
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin
        @require SymmetryReduceBZ="49a35663-c880-4242-bebb-1ec8c0fa8046" begin
            include("Visualization.jl")
            using .Visualization
            export Visualization
            export plotBZ
        end
    end
end

end
