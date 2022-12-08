
include("sodiumgreen.jl")
using .DFTGreen
using .DFTGreen.BrillouinZoneMeshes

using .BrillouinZoneMeshes.LinearAlgebra
using .BrillouinZoneMeshes.StaticArrays

using MCIntegration
using Random, Printf, BenchmarkTools, InteractiveUtils, Parameters

const Steps = 1e4

const scfres = DFTGreen.load_scfres("./run/sodium.jld2")
const gi = DFTGreen.GreenInterpolator(scfres)

const NGV = length(gi.gvectors)
# const NGV = 100
const τ = 0.0

function integrand(var, config)
    K, GV = var[1], var[2]
    # @assert idx == 1 "$(idx) is not a valid integrand"
    fk = SVector{3,Float64}(K[1], K[2], K[3])
    k = lattice_vector(gi.rbzmesh.mesh) * fk
    n1, n2 = GV[1], GV[2]
    gv1, gv2 = gi.gvectors[n1], gi.gvectors[n2]
    result = DFTGreen.green(gi, gv1, gv2, k, τ)
    return (real(result), imag(result))
    # return 1.0, 1.0
end

function measure(var, obs, weight, config)
    obs[1] += weight[1]
    obs[2] += weight[2]
end

function run(steps)

    K = MCIntegration.Continuous(-0.5, 0.5; alpha=2.0, adapt=true)
    Ext = MCIntegration.Discrete(1, NGV; adapt=true) # external variable is specified

    dof = [[3, 2], [3, 2]] # degrees of freedom of the normalization diagram and the bubble
    obs = zeros(Float64, 2)

    # config = MCIntegration.Configuration(var=(T, K, Ext), dof=dof, obs=obs, para=para)
    result = MCIntegration.integrate(integrand; measure=measure,
        var=(K, Ext), dof=dof, obs=obs, solver=:vegas,
        neval=steps, print=0, block=16)

    if isnothing(result) == false
        avg, std = result.mean, result.stdev

        println("avg=$avg, std=$std")
        # println(MCIntegration.summary(result))
        # i = 1
        # println(result.config.var[i].histogram)
        # println(sum(result.config.var[i].histogram))
        # println(result.config.var[i].accumulation)
        # println(result.config.var[i].distribution)
    end
end

run(Steps)
