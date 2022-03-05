module GridTree

using ..AbstractTrees
using ..StaticArrays
using ..BaseMesh

export GridNode, TreeGrid, uniformtreegrid, treegridfromdensity

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

function _calc_point(depth, pos, latvec)
    DIM = length(pos)
    ratio = pos ./ 2^depth .- 0.5
    origin = zeros(size(ratio))

    for i in 1:DIM
        origin .+= ratio .* latvec[i, :]
    end

    return origin
end

function _calc_origin(node::GridNode{DIM}, latvec) where {DIM}
    return _calc_point(node.depth, node.pos, latvec)
end

function _calc_cornerpoints(depth, pos, latvec)
    DIM = length(pos)
    points = []
    for i in 0:2^DIM-1
        ii = digits(i, base = 2, pad = DIM)
        push!(points, _calc_point(depth, pos .+ ii, latvec))
    end
    return points
end

function uniformtreegrid(isfine, latvec; maxdepth = 10, DIM = 2, N = 2) 
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

Base.length(tg::TreeGrid{DIM, SG}) where {DIM, SG} = length(tg.subgrids)
Base.size(tg::TreeGrid{DIM, SG}) where {DIM, SG} = length(tg) * size(tg.subgrids[1])
# index and iterator
function Base.getindex(tg::TreeGrid{DIM, SG}, i) where {DIM, SG}
    sgsize = size(tg.subgrids[1])

    tgi, sgi = (i - 1) ÷ sgsize + 1, (i - 1) % sgsize + 1

    return getindex(tg.subgrids[tgi], sgi)
end
Base.firstindex(tg::TreeGrid) = 1
Base.lastindex(tg::TreeGrid) = size(tg)

Base.iterate(tg::TreeGrid) = (tg[1],1)
Base.iterate(tg::TreeGrid, state) = (state>=size(tg)) ? nothing : (tg[state+1],state+1)

function densityisfine(density, latvec, depth, pos, rtol)
    cornerpoints = _calc_cornerpoints(depth, pos, latvec)
    cpval = [density(p) for p in cornerpoints]
    return abs(sum(cpval) / length(cpval) - density(_calc_point(depth, pos .+ 0.5, latvec))) < rtol
end

function treegridfromdensity(density, latvec; rtol = 1e-4, maxdepth = 10, DIM = 2, N = 2)
    isfine(depth, pos) = densityisfine(density, latvec, depth, pos, rtol)
    return uniformtreegrid(isfine, latvec; maxdepth = maxdepth, DIM = DIM, N = N)
end

end

