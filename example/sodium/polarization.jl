
include("sodiumgreen.jl")
using .DFTGreen
using .DFTGreen.BrillouinZoneMeshes

using .BrillouinZoneMeshes.LinearAlgebra
using .BrillouinZoneMeshes.StaticArrays

using MCIntegration
using Random, Printf, BenchmarkTools, InteractiveUtils, Parameters

const Steps = 1e6

const scfres = DFTGreen.load_scfres("./run/sodium.jld2")
const gi = DFTGreen.GreenInterpolator(scfres)

const NGV = length(gi.gvectors)
# const NGV = 100
const ω = 0.0
const gn1, gn2 = [0, 0, 0], [0, 0, 0]
const q = SVector{3,Float64}([0, 0, 0])

function integrand(var, config)
    T, K, GV = var[1], var[2], var[3]
    # @assert idx == 1 "$(idx) is not a valid integrand"
    τ = T[1]
    fk = SVector{3,Float64}(K[1], K[2], K[3])
    k2 = lattice_vector(gi.rbzmesh.mesh) * fk
    k1 = k2 .+ q
    gi1, gi2 = GV[1], GV[2]
    gm1, gm2 = gi.gvectors[gi1], gi.gvectors[gi2]
    result = DFTGreen.greenτ(gi, gn1 .+ gm1, gm1, k1, τ) * DFTGreen.greenτ(gi, gm1, gn2 .+ gm2, k2, gi.beta - τ)
    return result * exp(-im * ω * τ)
    # return 1.0, 1.0
end

function measure(var, obs, weight, config)
    obs[1] += weight[1]
end

function run(steps)

    T = MCIntegration.Continuous(0.0, gi.beta; alpha=2.0, adapt=true)
    K = MCIntegration.Continuous(-0.5, 0.5; alpha=2.0, adapt=true)
    Ext = MCIntegration.Discrete(1, NGV; adapt=true) # external variable is specified

    dof = [[1, 3, 2],] # degrees of freedom of the normalization diagram and the bubble
    obs = zeros(ComplexF64, 1)

    # config = MCIntegration.Configuration(var=(T, K, Ext), dof=dof, obs=obs, para=para)
    result = MCIntegration.integrate(integrand; measure=measure,
        var=(T, K, Ext), dof=dof, obs=obs, solver=:vegas,
        neval=steps, print=0, block=16, parallel=:thread, type=ComplexF64)

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
