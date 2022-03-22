using SpaceGrid
using SpaceGrid.AbstractTrees, SpaceGrid.GridTree, SpaceGrid.BaseMesh
using Plots
using LinearAlgebra

include("graphene.jl")

function density(k; band = 1)
    T = 0.01
    μ = -1.0

    dim = la.dim
    n = la.numatoms
    ϵ = dispersion(dim, n, ham, k)[band] - μ

    # return 1 / (exp((ϵ) / T) + 1.0)
    return (π*T)^2/((π*T)^2 + ϵ^2)
    # return (exp((ϵ) / T) / T)/(exp((ϵ) / T) + 1.0)^2
end

latvec = [2 0; 1 sqrt(3)] .* (2π)
# println(SpaceGrid.GridTree._calc_subpoints(3, [3,3], latvec, 2))
tg = treegridfromdensity(k->density(k), latvec; rtol = 1/2^14, maxdepth = 7, mindepth = 1, N = 2)

X, Y = zeros(Float64, size(tg)), zeros(Float64, size(tg))
for (pi, p) in enumerate(tg)
    X[pi], Y[pi] = p[1], p[2]
    # println(p)
end

println("size:$(size(tg)),\t length:$(length(tg)),\t efficiency:$(efficiency(tg))")

k = 4π
# p = plot(legend = false, size = (1024, 1024), xlim = (-1, 1), ylim = (-1, 1))
p = plot(legend = false, size = (1024, 1024), xlim = (-k, k), ylim = (-k, k))
scatter!(p, X, Y, marker = :cross, markersize = 2)
savefig(p, "run/tg.pdf")
