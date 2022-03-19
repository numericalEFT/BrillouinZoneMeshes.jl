using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh
using Plots
using LinearAlgebra

include("graphene.jl")

function density(k; band = 1)
    T = 0.001
    μ = -1.2

    dim = la.dim
    n = la.numatoms
    ϵ = dispersion(dim, n, ham, k)[band] - μ

    return 1 / (exp((ϵ) / T) + 1.0)
end

latvec = [2 0; 1 sqrt(3)] .* (2π)
tg = treegridfromdensity(k->density(k), latvec; rtol = 1e-3, maxdepth = 7, mindepth = 0, N = 2)

X, Y = zeros(Float64, size(tg)), zeros(Float64, size(tg))
for (pi, p) in enumerate(tg)
    X[pi], Y[pi] = p[1], p[2]
    # println(p)
end

println("size:$(size(tg))")
println("length:$(length(tg))")
println("efficiency:$(efficiency(tg))")

k = 4π
# p = plot(legend = false, size = (1024, 1024), xlim = (-1, 1), ylim = (-1, 1))
p = plot(legend = false, size = (1024, 1024), xlim = (-k, k), ylim = (-k, k))
scatter!(p, X, Y, marker = :cross, markersize = 2)
savefig(p, "run/tg.pdf")
