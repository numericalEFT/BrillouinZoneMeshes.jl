module DFTGreen

using DFTK
using Unitful
using UnitfulAtomic
using Interpolations
using BrillouinZoneMeshes
using JLD2

using BrillouinZoneMeshes.StaticArrays

include("interpolation.jl")
using .FFTInterp

export fourier_interpolate
export calc_scf, greenfromscf

struct GVectors{DIM} <: AbstractMeshes.AbstractMesh{Int,DIM}
    gmin::SVector{DIM,Int}
    # gmax::SVector{DIM,Int} gmax = gmin .+ size .- 1
    size::NTuple{DIM,Int}
end

function GVectors(gmin, gmax)
    DIM = length(gmin)
    @assert DIM == length(gmax)
    size = tuple((gmax .- gmin .+ 1)...)
    return GVectors{DIM}(gmin, size)
end

Base.getindex(gvs::GVectors, inds...) = gvs.gmin .+ inds .- 1
Base.getindex(gvs::GVectors, I::Int) = gvs[AbstractMeshes._ind2inds(size(gvs), I)...]

locate(gvs::GVectors, gv) = AbstractMeshes._inds2ind(size(gvs), (gv .- gvs.gmin .+ 1))
# compute scf, default sodium
function calc_scf(;
    lattice=(0.429u"nm") / 2 * [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]],
    pseudo_file=joinpath(@__DIR__, "Na_PBE.upf"),
    atoms=[ElementPsp(:Na, psp=load_psp(pseudo_file)),],
    positions=[[0.0, 0.0, 0.0],],
    smearing=DFTK.Smearing.Gaussian(),
    temperature=0.02u"Ry",
    Ecut=80.0u"Ry",
    kgrid=[16, 16, 16],
    tol=1e-8)

    model = model_PBE(lattice, atoms, positions; temperature, smearing)
    basis = PlaneWaveBasis(model; Ecut, kgrid=kgrid)
    scfres = self_consistent_field(basis, tol=tol)
    return scfres
end

function collect_gvector(scfres)
    gvectors = Vector{eltype(scfres.basis.G_vectors)}([])
    for (ik, k) in enumerate(scfres.basis.kpoints)
        for (ig, gv) in enumerate(k.G_vectors)
            if gv ∉ gvectors
                push!(gvectors, gv)
            end
        end
    end
    return gvectors
end

function gvector_minmax(gvectors)
    DIM = 3
    gmin = zeros(Int, DIM)
    gmax = zeros(Int, DIM)
    for (ig, gv) in enumerate(gvectors)
        for j in 1:DIM
            if gv[j] < gmin[j]
                gmin[j] = gv[j]
            end
            if gv[j] > gmax[j]
                gmax[j] = gv[j]
            end
        end
    end
    return gmin, gmax
end

function greenfromscf(scfres, m, n, kidx, τ, β)
    result = 0.0

    for band in 1:scfres.n_bands_converge
        ε = scfres.eigenvalues[kidx][band] - scfres.εF
        result += scfres.ψ[kidx][m, band] * conj(scfres.ψ[kidx][n, band]) * exp(-ε * τ) / (1 + exp(-ε * β))
    end
    return result
end

struct GreenInterpolator{DI}
    # store interpolators for green's function and provide interpolation
    # interpolators are for k-points, thus for each G-vector and band an interpolator is needed
    # store ψ and ε, green's function is computed
    beta::Float64
    n_bands::Int
    dispersions::Vector{DI}
    ψ
    kgrid
    gvectors::GVectors{3}
    rbzmesh
end

function GreenInterpolator(scfres)
    # extract parameters
    beta = scfres.εF / scfres.basis.model.temperature
    n_bands = scfres.n_bands_converge
    kgrid = scfres.basis.kgrid
    lattice = scfres.basis.model.lattice
    atoms = ones(Int, length(scfres.basis.model.atoms))
    positions = scfres.basis.model.positions

    # generate bzmesh and map
    cell = Cells.Cell(; lattice=Matrix(austrip.(lattice)), atoms=atoms, positions=positions)
    bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid, shift=[false, false, false])
    meshmap = MeshMaps.MeshMap(bzmesh)
    rbzmesh = MeshMaps.ReducedBZMesh(bzmesh, meshmap)

    kx = LinRange(-0.5, 0.5 - 1.0 / kgrid[1], kgrid[1])
    ky = LinRange(-0.5, 0.5 - 1.0 / kgrid[2], kgrid[2])
    kz = LinRange(-0.5, 0.5 - 1.0 / kgrid[3], kgrid[3])


    gvfromscf = collect_gvector(scfres)
    gmin, gmax = gvector_minmax(gvfromscf)
    gvectors = GVectors(gmin, gmax)

    # collect data

    ψ = zeros(ComplexF64, (length(gvectors), scfres.n_bands_converge, length(rbzmesh)))
    ϵ = zeros(Float64, (scfres.n_bands_converge, length(rbzmesh)))
    for bandidx in 1:scfres.n_bands_converge
        for ki in 1:length(scfres.basis.kpoints)
            k = scfres.basis.kpoints[ki]
            idx = AbstractMeshes.locate(rbzmesh, lattice_vector(rbzmesh.mesh) * k.coordinate)
            ϵ[bandidx, idx] = scfres.eigenvalues[ki][bandidx] - scfres.εF
            for gi in 1:length(k.G_vectors)
                g = k.G_vectors[gi]
                gidx = locate(gvectors, g)
                @assert gvectors[gidx] ≈ g
                ψ[gidx, bandidx, idx] = scfres.ψ[ki][gi, bandidx]
            end
        end
    end

    # generate interpolator

    dispersions = []

    for bandidx in 1:scfres.n_bands_converge
        # dispersion
        bandarray = zeros(kgrid...)
        for kx in 1:kgrid[1]
            for ky in 1:kgrid[2]
                for kz in 1:kgrid[3]
                    # idx = AbstractMeshes._inds2ind(tuple(kgrid...), [kx, ky, kz])
                    idx = AbstractMeshes.locate(rbzmesh, rbzmesh[kx, ky, kz])
                    bandarray[kx, ky, kz] = ϵ[bandidx, AbstractMeshes.locate(meshmap, idx)]
                end
            end
        end
        itp = interpolate(bandarray, BSpline(Cubic(Line(OnGrid()))))
        sitp = scale(itp, kx, ky, kz)
        push!(dispersions, sitp)
    end

    return GreenInterpolator{typeof(dispersions[1])}(beta, n_bands, dispersions,
        ψ, kgrid, gvectors, rbzmesh)
end

function green(gi::GreenInterpolator{DI}, gv1, gv2, k, τ) where {DI}
    gi1, gi2 = locate(gi.gvectors, gv1), locate(gi.gvectors, gv2)
    result = ComplexF64(0.0)

    for band in 1:gi.n_bands
        data1 = view(gi.ψ, gi1, band, :)
        ψ1 = AbstractMeshes.interp(data1, gi.rbzmesh, k)
        data2 = view(gi.ψ, gi2, band, :)
        ψ2 = AbstractMeshes.interp(data2, gi.rbzmesh, k)
        sitp = gi.dispersions[band]
        ε = sitp(k...)
        # println("ψ1=$ψ1, ψ2=$ψ2, ε=$ε")
        result += ψ1 * conj(ψ2) * exp(-ε * τ) / (1 + exp(-ε * gi.beta))
    end
    return result
end


end



if abspath(PROGRAM_FILE) == @__FILE__

end