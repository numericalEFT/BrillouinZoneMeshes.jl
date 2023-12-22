
# test barycheb integrate on an arbitrary region

using Test
using Random

using BrillouinZoneMeshes

using BrillouinZoneMeshes.BaseMesh
using BrillouinZoneMeshes.LinearAlgebra
using BrillouinZoneMeshes.CompositeGrids.BaryChebTools
using BrillouinZoneMeshes.StaticArrays

function interpweight(n, xs, wc, xc, DIM)
    inds = CartesianIndices(NTuple{DIM,Int}(ones(Int, DIM) .* n))
    haseq = false
    eqinds = zeros(Int, DIM)
    for i in 1:DIM
        for j in 1:n
            if xs[i] == xc[j]
                eqinds[i] = j
                haseq = true
            end
        end
    end

    result = zeros(Float64, size(inds))

    if haseq
        newxs = [xs[i] for i in 1:DIM if eqinds[i] == 0]
        newDIM = length(newxs)
        if newDIM == 0
            result[CartesianIndex(eqinds...)] = 1.0
            return result
        else
            newresult = view(result, [(i == 0) ? (1:n) : (i) for i in eqinds]...)
            newresult .= _interpweight_noneq(n, newxs, wc, xc, newDIM)
            return result
        end
    else
        return _interpweight_noneq(n, xs, wc, xc, DIM)
    end
end

function _interpweight_noneq(n, xs, wc, xc, DIM)
    # deal with the case when there's no xs[i] = xc[j]
    inds = CartesianIndices(NTuple{DIM,Int}(ones(Int, DIM) .* n))
    den = 0.0
    result = zeros(Float64, size(inds))
    for (indi, ind) in enumerate(inds)
        q = 1.0
        for i in 1:DIM
            q *= wc[ind[i]] / (xs[i] - xc[ind[i]])
        end
        result[indi] += q
        den += q
    end
    return result ./ den
end

function interpweightND(xgrid::BaryCheb1D{N}, xs) where {N}
    return interpweight(N, xs, xgrid.w, xgrid.x, length(xs))
end

function interpweight(mesh::ChebMesh{T,DIM,N}, x) where {T,DIM,N}
    displacement = SVector{DIM,Float64}(x) - mesh.origin
    xs = AbstractMeshes.cart_to_frac(mesh, displacement) .* 2.0 .- 1.0
    return interpweightND(mesh.barycheb, xs)
end

@testset "BaryCheb integrate on arbitrary region" begin
    rng = MersenneTwister(1234)

    radius = 1.0

    function region(x)
        # define region
        # return 1.0 if in region, 0.0 if out

        # a circle
        if dot(x, x) < radius
            return 1.0
        else
            return 0.0
        end
    end

    function integrand(x)
        return exp(-dot(x, x))
    end

    N = 10
    origin = [-1, -1.0]
    lattice = [2 0; 0 2.0]
    cm = BaseMesh.ChebMesh(origin, lattice, 2, N)
    barycheb = cm.barycheb
    data = zeros(Float64, size(cm))
    for (i, p) in enumerate(cm)
        data[i] = integrand(p)
    end

    @testset "interp weight" begin
        for i in 1:10
            x, y = rand(rng) * 2 - 1, rand(rng) * 2 - 1
            wmat = interpweight(cm, [x, y])
            @test isapprox(sum(wmat .* data), integrand([x, y]), rtol=1e-3)
        end
    end

    @testset "integrate on region" begin
        # record weight 
        M = 200
        um = BaseMesh.UMesh{Float64,2}(cm.lattice, cm.inv_lattice, cm.cell_volume, cm.origin, (M, M), [0.0, 0.0])
        weight = zeros(Float64, size(cm))
        for (i, p) in enumerate(um)
            weight .+= region(p) .* interpweight(cm, p) .* AbstractMeshes.volume(um, i)
        end

        # compute
        result = sum(weight .* data)
        println("result=$(result), acc. = $((ℯ - 1) / ℯ * π)")
        @test isapprox(result, (ℯ - 1) / ℯ * π, rtol=1e-2)
        @test isapprox(sum(weight), π, rtol=1e-2)
    end



end