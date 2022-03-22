module GridTree

using ..AbstractTrees
using ..StaticArrays
using ..BaseMesh
using ..Statistics
using ..LinearAlgebra
export GridNode, TreeGrid, uniformtreegrid, treegridfromdensity, efficiency

struct GridNode{DIM}
    depth::Int
    pos::SVector{DIM, Int64}
    children::Vector{GridNode{DIM}}
end

function GridNode{DIM}(
    isfine;
    depth = 0, pos = SVector{DIM, Int64}(zeros(Int64, DIM)), maxdepth = 10, mindepth = 0) where {DIM}
    if (isfine(depth, pos) && depth >= mindepth) || depth >= maxdepth
        return GridNode{DIM}(depth, pos, Vector{GridNode{DIM}}([]))
    else
        children = Vector{GridNode{DIM}}([])
        for i in 0:2^DIM-1
            bin = digits(i, base = 2, pad = DIM) |> reverse
            childdepth = depth + 1
            childpos = pos .* 2 .+ bin
            push!(children, GridNode{DIM}(isfine; depth = childdepth, pos = childpos, maxdepth = maxdepth, mindepth = mindepth))
        end
        return GridNode{DIM}(depth, pos, children)
    end
end

AbstractTrees.children(node::GridNode) = node.children

function efficiency(root::GridNode{DIM}) where {DIM}
    np, depth = 0, 0
    for node in PostOrderDFS(root)
        if node.depth > depth
            depth = node.depth
        end
        if isempty(node.children)
            np = np + 1
        end
    end
    return np / 2^(depth*DIM)
end

struct TreeGrid{DIM, SG}
    root::GridNode{DIM}
    latvec::SMatrix{DIM, DIM, Float64}
    subgrids::Vector{SG}
end

efficiency(tg::TreeGrid) = efficiency(tg.root)

function _calc_area(latvec)
    return abs(det(latvec))
end

function _calc_point(depth, pos, latvec)
    DIM = length(pos)
    ratio = pos ./ 2^depth .- 0.5
    origin = zeros(size(ratio))

    for i in 1:DIM
        origin .+= ratio[i] .* latvec[i, :]
    end

    return origin
end

function _calc_origin(node::GridNode{DIM}, latvec) where {DIM}
    return _calc_point(node.depth, node.pos, latvec)
end

function _calc_subpoints(depth, pos, latvec, N)
    DIM = length(pos)
    points = []
    for i in 0:N^DIM-1
        ii = digits(i, base = N, pad = DIM)
        push!(points, _calc_point(depth, pos .+ ii / (N-1), latvec))
    end
    return points
end

function _calc_cornerpoints(depth, pos, latvec)
    # DIM = length(pos)
    # points = []
    # for i in 0:2^DIM-1
    #     ii = digits(i, base = 2, pad = DIM)
    #     push!(points, _calc_point(depth, pos .+ ii, latvec))
    # end
    # return points
    return _calc_subpoints(depth, pos, latvec, 2)
end

function uniformtreegrid(isfine, latvec; maxdepth = 10, mindepth = 0, DIM = 2, N = 2)
    root = GridNode{DIM}(isfine; maxdepth = maxdepth, mindepth = mindepth)
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

function densityisfine(density, latvec, depth, pos, rtol, DIM; N = 3)
    # compare results from subgrid of N+1 and N+3
    # area = _calc_area(latvec) / 2^(depth*DIM)
    area = 1.0 / 2^(depth*DIM)
    cornerpoints1 = _calc_subpoints(depth, pos, latvec, N)
    cornerpoints2 = _calc_subpoints(depth, pos, latvec, 16)
    val1 = [density(p) for p in cornerpoints1]
    val2 = [density(p) for p in cornerpoints2]
    # return abs(sum(val1) / length(val1) - sum(val2) / length(val2)) * area < rtol
    # return abs(sum(val2) / length(val2)) * area < rtol
    # println("max:$(abs(maximum(val2)))")
    # println("area:$(area)")
    # println("rtol:$(rtol)")
    return abs(maximum(val2)) * area < rtol
    # return std(val2) * area < rtol
end

function treegridfromdensity(density, latvec; rtol = 1e-4, maxdepth = 10, mindepth = 0, DIM = 2, N = 2)
    isfine(depth, pos) = densityisfine(density, latvec, depth, pos, rtol, DIM; N = N)
    return uniformtreegrid(isfine, latvec; maxdepth = maxdepth, mindepth = mindepth, DIM = DIM, N = N)
end

end

