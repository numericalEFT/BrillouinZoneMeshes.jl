module DFTGreen

using DFTK
using Unitful
using UnitfulAtomic
using Interpolations
using BrillouinZoneMeshes

include("interpolation.jl")
using .FFTInterp

export fourier_interpolate
export calc_scf

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

end


if abspath(PROGRAM_FILE) == @__FILE__
    # using .DFTGreen

    using .DFTGreen.DFTK
    using .DFTGreen.Unitful
    using .DFTGreen.UnitfulAtomic
    using .DFTGreen.Interpolations
    using .DFTGreen.BrillouinZoneMeshes

    using Plots
    # using PlotlyJS
    using Printf

    scfres = calc_scf()

    function greenfromscf(scfres, m, n, kidx, τ, β)
        result = 0.0

        for band in 1:scfres.n_bands_converge
            ε = scfres.eigenvalues[kidx][band] - scfres.εF
            result += scfres.ψ[kidx][m, band] * conj(scfres.ψ[kidx][n, band]) * exp(-ε * τ) / (1 + exp(-ε * β))
        end
        return result
    end

    # build BZ mesh and maps
    cell = Cells.Cell(; lattice=austrip.(lattice), atoms=[1,], positions)
    bzmesh = BZMeshes.Monkhorst_Pack(cell=cell, size=kgrid, shift=[false, false, false])
    meshmap = MeshMaps.MeshMap(bzmesh)

    println(length(meshmap.irreducible_indices))
    for idx in meshmap.irreducible_indices
        k = inv_lattice_vector(bzmesh) * bzmesh[idx]
        println(k)
    end

    # band structure with BZ mesh and maps
    band = Dict{Int,Float64}()
    bandix = 5 # the fifth band is the valence band
    for (ki, kpoints) in enumerate(basis.kpoints)
        idx = AbstractMeshes.locate(bzmesh, lattice_vector(bzmesh) * basis.kpoints[ki].coordinate)
        band[meshmap.map[idx]] = greenfromscf(scfres, 1, 1, ki, 0.0, 100.0)
    end

    # band fourier transform
    bandarray = zeros(kgrid...)

    for kx in 1:kgrid[1]
        for ky in 1:kgrid[2]
            for kz in 1:kgrid[3]
                idx = AbstractMeshes._inds2ind(tuple(kgrid...), [kx, ky, kz])
                # idx = AbstractMeshes.locate(bzmesh, bzmesh[kx, ky, kz])
                bandarray[kx, ky, kz] = band[meshmap.map[idx]]
            end
        end
    end

    kx = LinRange(-0.5, 0.5 - 1.0 / kgrid[1], kgrid[1])
    ky = LinRange(-0.5, 0.5 - 1.0 / kgrid[2], kgrid[2])
    kz = LinRange(-0.5, 0.5 - 1.0 / kgrid[3], kgrid[3])
    itp = interpolate(bandarray, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, kx, ky, kz)

    Kfine = 100
    Kfinegrid = collect(LinRange(kx[1], kx[end], Kfine))

    p = Plots.plot()
    kyz = [(1, 1), (4, 4), (8, 8), (12, 12), (16, 16)]
    # kyi, kzi = Int(kgrid[2] / 2), Int(kgrid[3] / 2)
    for (kyi, kzi) in kyz
        Plots.scatter!(p, kx, bandarray[:, kyi, kzi], label="band (ky=$kyi, kz = $kzi)", markershape=:circle)
        band_interpolate = [sitp(k, ky[kyi], kz[kzi]) for k in Kfinegrid]
        Plots.plot!(p, Kfinegrid, band_interpolate, label="interpolation")
    end
    display(p)

end