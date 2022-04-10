using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh
using LinearAlgebra
using Test

@testset "SpaceGrid.jl" begin
    if isempty(ARGS)
        include("tree.jl")
        include("basemesh.jl")
    else
        include(ARGS[1])
    end
end
