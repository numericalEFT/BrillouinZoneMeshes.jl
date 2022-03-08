using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh
using Plots
using LinearAlgebra

mesh = UniformMesh{2, 4}([0.0, 0.0], [0 1 ; 1 0])
show(mesh)

naiveisfine(depth, pos) = depth >= 2

# tree = GridNode{2}(naiveisfine)
# print_tree(tree)

# for node in PostOrderDFS(tree)
#     if isempty(node.children)
#         println(node.pos)
#     end
# end

function dispersion(k)
    me = 0.5
    μ = 1.0
    return norm(k)^2/(2me)-μ

    # t = 0.5
    # μ = 0.0
    # return sum([t * cos(π*p) for p in k]) + μ
 end

function density(K, latvec)
    DIM = length(K)

    T = 0.001

    kpoints = [K, ]
    for i in 1:DIM
        push!(kpoints, K .+ latvec[i, :])
        push!(kpoints, K .- latvec[i, :])
    end

    ϵ = minimum([dispersion(k) for k in kpoints])

    return 1 / (exp((ϵ) / T) + 1.0)
    # return 1 / ((π * T)^2 + ϵ^2)
end

latvec = [2 0; 1 sqrt(3)]
# latvec = [2 0; 0 2]
# tg = uniformtreegrid(naiveisfine, latvec; N = 3)
# tg = treegridfromdensity(density, latvec; rtol = 1e-4, maxdepth = 5)
tg = treegridfromdensity(k->density(k, latvec), latvec; rtol = 1e-1, maxdepth = 6)

X, Y = zeros(Float64, size(tg)), zeros(Float64, size(tg))
for (pi, p) in enumerate(tg)
    X[pi], Y[pi] = p[1], p[2]
    # println(p)
end

println("size:$(size(tg))")
println("length:$(length(tg))")
println("efficiency:$(efficiency(tg))")

# p = plot(legend = false, size = (1024, 1024), xlim = (-1, 1), ylim = (-1, 1))
p = plot(legend = false, size = (1024, 1024), xlim = (-1.5, 1.5), ylim = (-1.5, 1.5))
scatter!(p, X, Y, marker = :cross, markersize = 2)
savefig(p, "run/tg.pdf")
