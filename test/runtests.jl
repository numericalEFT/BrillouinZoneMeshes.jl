module _Test_BrillouinZoneMeshes

using BrillouinZoneMeshes
using BrillouinZoneMeshes.AbstractTrees
# using BrillouinZoneMeshes.TreeMeshes
# using BrillouinZoneMeshes.BaseMesh
# using BrillouinZoneMeshes.BaryCheb
# using BrillouinZoneMeshes.SymMaps

using LinearAlgebra, Random
using Test

include("testcase.jl")

@testset "BrillouinZoneMeshes.jl" begin

    if isempty(ARGS)
        # include("barycheb.jl")
        include("BaseMesh.jl")
        include("TreeMeshes.jl")
        include("mc.jl")
        include("PointSymmetry.jl")
        include("UniformMeshMap.jl")
        include("MeshMap.jl")
        include("Model.jl")
    else
        include(ARGS[1])
    end
end

end

