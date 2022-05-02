using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh, SpaceGrid.BaryCheb
using LinearAlgebra, Random
using Test

@testset "SpaceGrid.jl" begin
    if isempty(ARGS)
        include("tree.jl")
        include("basemesh.jl")
        include("barycheb.jl")
    else
        include(ARGS[1])
    end
end
