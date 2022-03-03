module SpaceGrid

using AbstractTrees
using StaticArrays

# Write your package code here.

include("basemesh.jl")
using .BaseMesh

include("tree.jl")
using .GridTree

end
