using BrillouinZoneMeshes
using BrillouinZoneMeshes.AbstractTrees, BrillouinZoneMeshes.TreeMeshes, BrillouinZoneMeshes.BaseMesh, BrillouinZoneMeshes.BaryCheb
using LinearAlgebra, Random
using Test

@testset "BrillouinZoneMeshes.jl" begin


    if isempty(ARGS)
        include("barycheb.jl")
        include("BaseMesh.jl")
        include("TreeMeshes.jl")
        include("mc.jl")
    else
        include(ARGS[1])
    end
end

