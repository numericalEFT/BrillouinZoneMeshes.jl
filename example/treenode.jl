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
    μ = 0.2
    return norm(k)^2/(2me)-μ

    # t = 0.5
    # μ = 0.0
    # return sum([t * cos(π*p) for p in k]) + μ
 end

function density(K, latvec)
    DIM = length(K)

    T = 0.01

    kpoints = [K, ]
    for i in 0:3^DIM-1
        ii = digits(i, base = 3, pad = DIM)
        Knew = copy(K)
        for j in 1:DIM
            Knew = Knew .+ (ii[j] - 1) .* latvec[j, :]
        end
        push!(kpoints, Knew)
    end

    Kmin = kpoints[findmin([norm(k) for k in kpoints])[2]]

    # ϵ = minimum([dispersion(k) for k in kpoints])
    ϵ = dispersion(Kmin)

    return 1 / (exp((ϵ) / T) + 1.0)
    # return (π * T)^2 / ((π * T)^2 + ϵ^2)
end

latvec = [2 0; 1 sqrt(3)]
# latvec = [2 0; 0 2]
# tg = uniformtreegrid(naiveisfine, latvec; N = 3)
# tg = treegridfromdensity(density, latvec; rtol = 1e-4, maxdepth = 5)
tg = treegridfromdensity(k->density(k, latvec), latvec; rtol = 1e-3, maxdepth = 7, mindepth = 2, N = 2)

X, Y = zeros(Float64, size(tg)), zeros(Float64, size(tg))
for (pi, p) in enumerate(tg)
    X[pi], Y[pi] = p[1], p[2]
    # println(p)
end

println("size:$(size(tg))")
println("length:$(length(tg))")
println("efficiency:$(efficiency(tg))")

# p = plot(legend = false, size = (1024, 1024), xlim = (-1, 1), ylim = (-1, 1))
p = plot(legend = false, size = (1024, 1024), xlim = (-2.0, 2.0), ylim = (-2.0, 2.0))
scatter!(p, X, Y, marker = :cross, markersize = 2)
savefig(p, "run/tg.pdf")
