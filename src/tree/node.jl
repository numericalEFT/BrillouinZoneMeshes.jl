module GridTree

using ..AbstractTrees
using ..StaticArrays

export GridNode

struct GridNode{DIM}
    depth::Int
    pos::SVector{DIM, Int64}
    children::Vector{GridNode{DIM}}
end

AbstractTrees.children(node::GridNode) = node.children

function GridNode{DIM}(isfine; depth = 0, pos = SVector{DIM, Int64}(zeros(Int64, DIM)), maxdepth = 10) where {DIM}
    if isfine(depth, pos) || depth == maxdepth
        return GridNode{DIM}(depth, pos, [])
    else
        children = Vector{GridNode{DIM}}([])
        for i in 0:2^DIM-1
            bin = digits(i, base = 2, pad = DIM) |> reverse
            childdepth = depth + 1
            childpos = pos .* 2 .+ bin
            push!(children, GridNode{DIM}(isfine; depth = childdepth, pos = childpos, maxdepth = 10))
        end
        return GridNode{DIM}(depth, pos, children)
    end
end

end

