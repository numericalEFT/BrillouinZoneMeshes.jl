
include("sodiumgreen.jl")
using .DFTGreen
using .DFTGreen.BrillouinZoneMeshes

using .BrillouinZoneMeshes.LinearAlgebra
using .BrillouinZoneMeshes.StaticArrays

using MCIntegration
using Random, Printf, BenchmarkTools, InteractiveUtils, Parameters
using Plots

const Steps = 1e6

const scfres = DFTGreen.load_scfres("./run/sodium.jld2")
const gi = DFTGreen.GreenInterpolator(scfres)

const NGV = length(gi.gvectors)
# const NGV = 100
const δ = 1e-4

const NExt = 20
const ωmin, ωmax = -0.07, 0.15
const ωlist = LinRange(ωmin, ωmax, NExt)

function integrand(var, config)
    Ext, K, GV = var[1], var[2], var[3]
    # @assert idx == 1 "$(idx) is not a valid integrand"
    ω = ωlist[Ext[1]]
    fk = SVector{3,Float64}(K[1], K[2], K[3])
    k = lattice_vector(gi.rbzmesh.mesh) * fk
    n1, n2 = GV[1], GV[2]
    gv1, gv2 = gi.gvectors[n1], gi.gvectors[n2]
    result = DFTGreen.greenω(gi, gv1, gv2, k, ω; δ=δ)
    return (real(result), imag(result))
    # return 1.0, 1.0
end

function measure(var, obs, weight, config)
    n = var[1][1]
    obs[1][n] += weight[1]
    obs[2][n] += weight[2]
end

function run(steps)

    Ext = MCIntegration.Discrete(1, NExt; adapt=true)
    K = MCIntegration.Continuous(-0.5, 0.5; alpha=2.0, adapt=true)
    GV = MCIntegration.Discrete(1, NGV; adapt=true)

    dof = [[1, 3, 2], [1, 3, 2]] # degrees of freedom of the normalization diagram and the bubble
    obs = [zeros(Float64, NExt), zeros(Float64, NExt)]

    # config = MCIntegration.Configuration(var=(T, K, Ext), dof=dof, obs=obs, para=para)
    result = MCIntegration.integrate(integrand; measure=measure,
        var=(Ext, K, GV), dof=dof, obs=obs, solver=:vegas,
        neval=steps, print=0, block=16, parallel=:thread)

    if isnothing(result) == false
        avg, std = result.mean, result.stdev

        println("avg=$avg, std=$std")
        # println(MCIntegration.summary(result))
        # i = 1
        # println(result.config.var[i].histogram)
        # println(sum(result.config.var[i].histogram))
        # println(result.config.var[i].accumulation)
        # println(result.config.var[i].distribution)
        plt = plot(ωlist ./ scfres.εF, -avg[2], yerror=std[2], xlabel="ω", label="circle")#, xlims=[ωmin, ωmax])
        display(plt)
        readline()
        savefig(plt, "run/dos.png")
        return avg, std
    end
end

run(Steps)
