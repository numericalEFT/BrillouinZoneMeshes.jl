using BZMeshes
using BZMeshes.AbstractTrees, BZMeshes.GridTree, BZMeshes.BaseMesh, BZMeshes.BaryCheb
using LinearAlgebra, Random
using Test

@testset "BZMeshes.jl" begin


    if isempty(ARGS)
        include("tree.jl")
        include("basemesh.jl")
        include("barycheb.jl")
        include("mc.jl")
    else
        include(ARGS[1])
    end
end

