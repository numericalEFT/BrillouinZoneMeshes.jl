using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh
using Plots
using LinearAlgebra

mesh = UniformMesh{2, 4}([0.0, 0.0], [0 1 ; 1 0])
show(mesh)

naiveisfine(depth, pos) = depth >= 2

tree = GridNode{2}(naiveisfine)
print_tree(tree)

for node in PostOrderDFS(tree)
    if isempty(node.children)
        println(node.pos)
    end
end

function density(K)
    me = 0.5
    T = 0.01
    μ = 1.0

    k = norm(K)
    ϵ = k^2 / (2me)

    return 1 / (exp((ϵ - μ) / T) + 1.0)
end

# tg = uniformtreegrid(naiveisfine, [0 1; 1 0])
tg = treegridfromdensity(density, [0 2; 2 0]; rtol = 5e-2, maxdepth = 5)

X, Y = zeros(Float64, size(tg)), zeros(Float64, size(tg))
for (pi, p) in enumerate(tg)
    X[pi], Y[pi] = p[1], p[2]
    # println(p)
end

println("size:$(size(tg))")
println("length:$(length(tg))")

p = plot(legend = false, size = (1024, 1024))
scatter!(p, X, Y, marker = :cross, markersize = 2)
savefig(p, "run/tg.pdf")
