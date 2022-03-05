module GridTree

using ..AbstractTrees
using ..StaticArrays
using ..BaseMesh

export GridNode, TreeGrid, uniformtreegrid

struct GridNode{DIM}
    depth::Int
    pos::SVector{DIM, Int64}
    children::Vector{GridNode{DIM}}
end
function GridNode{DIM}(isfine; depth = 0, pos = SVector{DIM, Int64}(zeros(Int64, DIM)), maxdepth = 10) where {DIM}
    if isfine(depth, pos) || depth == maxdepth
        return GridNode{DIM}(depth, pos, Vector{GridNode{DIM}}([]))
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


AbstractTrees.children(node::GridNode) = node.children

struct TreeGrid{DIM, SG}
    root::GridNode{DIM}
    latvec::SMatrix{DIM, DIM, Float64}
    subgrids::Vector{SG}
end

function _calc_origin(node::GridNode{DIM}, latvec) where {DIM}
    ratio = node.pos ./ 2^node.depth
    origin = zeros(size(ratio))

    for i in 1:DIM
        origin .+= ratio .* latvec[i, :]
    end

    return origin
end

function uniformtreegrid(isfine, latvec; maxdepth = 10, DIM=2, N=2) 
    root = GridNode{DIM}(isfine; maxdepth = maxdepth)
    subgrids = Vector{UniformMesh{DIM, N}}([])
    for node in PostOrderDFS(root)
        if isempty(node.children)
            depth = node.depth
            origin = _calc_origin(node, latvec )
            mesh = UniformMesh{DIM, N}(origin, latvec ./ 2^depth)
            push!(subgrids, mesh)
        end
    end
    return TreeGrid{DIM, UniformMesh{DIM, N}}(root, latvec, subgrids)
end

end

