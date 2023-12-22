
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

### extract parameters
# const beta = 10.0
const DIM = 3
const beta = gi.beta
const EF = scfres.εF
const kF = 0.4795 # an approx val, used for sampling
const kmax = norm(DFTGreen.lattice_vector(gi.rbzmesh.mesh) * [1, 1, 1]) * 1.1 # ∀k ∈ 1st BZ, |kmax|>|k|

const NGV = length(gi.gvectors)
# const NGX = maximum(size(gi.gvectors))
const NGX = 21
const NMIN = 11

# const NGV = 100
const ω = 0.0
const gn1, gn2 = SVector{3,Int}(0, 0, 0), SVector{3,Int}(0, 0, 0)
const q = SVector{3,Float64}([0, 0, 0])

println("beta=$(beta)")
println("NGV=$(NGV)")
println("lattice = $(lattice_vector(gi.rbzmesh.mesh))")

function is_in_bz(bzmesh, k)
    fk = DFTGreen.inv_lattice_vector(bzmesh) * k
    return all((-0.5 <= fki <= 0.5) for fki in fk)
end

function integrand(var, config)
    # factor = cell_volume(gi.rbzmesh.mesh)
    factor = 1.0
    T, K, GV = var[1], var[2], var[3]
    # @assert idx == 1 "$(idx) is not a valid integrand"
    τ = T[1]

    # fk = SVector{3,Float64}(K[1], K[2], K[3])
    # k2 = lattice_vector(gi.rbzmesh.mesh) * fk
    # k1 = k2 .+ q
    k2 = K[1]
    if !is_in_bz(gi.rbzmesh.mesh, k2)
        return ComplexF64(0.0)
    end
    k1 = k2 .+ q

    gm1 = SVector{3,Int}(GV[1] - NMIN, GV[2] - NMIN, GV[3] - NMIN)
    gm2 = SVector{3,Int}(GV[4] - NMIN, GV[5] - NMIN, GV[6] - NMIN)
    # gi1, gi2 = GV[1], GV[2]
    # gm1, gm2 = gi.gvectors[gi1], gi.gvectors[gi2]
    result = -DFTGreen.greenτ(gi, gn1 .+ gm1, gm2, k1, τ; beta=beta) * DFTGreen.greenτ(gi, gn2 .+ gm2, gm1, k2, beta - τ; beta=beta)
    return result * exp(-im * ω * τ) * factor
    # return 1.0, 1.0
end

function integrand2(var, config)
    factor = cell_volume(gi.rbzmesh.mesh) / (2π)^3
    # UEG, greenfree, gn1=gn2=[0 0 0], gi1==gi2, q=0
    T, K, GV = var[1], var[2], var[3]
    # @assert idx == 1 "$(idx) is not a valid integrand"
    τ = T[1]
    fk = SVector{3,Float64}(K[1], K[2], K[3])
    k2 = lattice_vector(gi.rbzmesh.mesh) * fk
    k1 = k2 .+ q
    # gi1 = GV[1]
    gv1 = [GV[1] - NMIN, GV[2] - NMIN, GV[3] - NMIN]
    gv2 = [GV[4] - NMIN, GV[5] - NMIN, GV[6] - NMIN]
    # gm1 = gi.gvectors[gi1]
    # result = -DFTGreen.greenfreeτ(gi, gi1, gi1, k1, τ; beta=beta) * DFTGreen.greenfreeτ(gi, gi1, gi1, k2, beta - τ; beta=beta)
    result = -DFTGreen.greenfreeτ(gi, gv1, gv2, k1, τ; beta=beta) * DFTGreen.greenfreeτ(gi, gv1, gv2, k2, beta - τ; beta=beta)
    return result * exp(-im * ω * τ) * factor
end

function measure(var, obs, weight, config)
    obs[1] += weight[1]
end

function run(steps)

    T = MCIntegration.Continuous(0.0, beta; alpha=2.0, adapt=true)
    # K = MCIntegration.Continuous(-0.5, 0.5; alpha=2.0, adapt=true)
    K = MCIntegration.FermiK(3, kF, 0.2kF, kmax)
    # Ext = MCIntegration.Discrete(1, NGV; adapt=true) # external variable is specified
    Ext = MCIntegration.Discrete(1, NGX; adapt=true)

    dof = [[1, 1, 6],] # degrees of freedom of the normalization diagram and the bubble
    # dof = [[1, 3, 2],] # degrees of freedom of the normalization diagram and the bubble
    obs = zeros(ComplexF64, 1)

    # config = MCIntegration.Configuration(var=(T, K, Ext), dof=dof, obs=obs, para=para)
    result = MCIntegration.integrate(integrand; measure=measure,
        var=(T, K, Ext), dof=dof, obs=obs, solver=:vegasmc,
        neval=steps / 10, print=0, block=16, parallel=:thread, type=ComplexF64)
    result = MCIntegration.integrate(integrand; measure=measure,
        config=result.config,
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
