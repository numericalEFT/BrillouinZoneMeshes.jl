using DFTK
using Unitful
using UnitfulAtomic
using Plots
using PlotlyJS
using Printf

# Setup a sodium lattice, see the following link for more information
# http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/bcc/bcc.php
a = 0.429u"nm"

# each column is a lattice vector
lattice = a / 2 * [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]]
pseudo_file = joinpath(@__DIR__, "Na_PBE.upf")
Na = ElementPsp(:Na, psp=load_psp(pseudo_file))
atoms = [Na,]
positions = [[0.0, 0.0, 0.0],]

smearing = DFTK.Smearing.Gaussian()
temperature = 0.02u"Ry"

model = model_PBE(lattice, atoms, positions; temperature, smearing)

Ecut = 80.0u"Ry"
kgrid = [8, 8, 8]

basis = PlaneWaveBasis(model; Ecut, kgrid=kgrid)

scfres = self_consistent_field(basis, tol=1e-8);

# plot
plot_dos(scfres, unit=u"eV", xlims=[-5, 20], ylims=[0, 40]) # in hartree unit
plot_bandstructure(scfres, kline_density=10, unit=u"eV", ρ=scfres.ρ)

# calculate band width
_, data = compute_bands(basis, [[0.0, 0.0, 0.0],], ρ=scfres.ρ);

bandwidth = uconvert.(u"eV", (data[1][5] - scfres.εF) * u"hartree")
println(bandwidth)

# analysis the wavefunction
# the 1s band is the 2s^2 orbital
# the 2-4 bands are the 2p^6 orbitals
# the >5 bands are the valence electrons

kidx = 1

for (gi, g) in enumerate(scfres.basis.kpoints[kidx].G_vectors)
    ψ=scfres.ψ[kidx]
    @printf("%32s    %16.8f    %16.8f    %16.8f\n", "$g", abs(ψ[gi, 1]), abs(ψ[gi, 2]), abs(ψ[gi, 3])) 
end

Gz = 0
bandidx = 1 

x, y, psi = [], [], []

for (gi, g) in enumerate(scfres.basis.kpoints[kidx].G_vectors)
    if g[3]==Gz
        push!(x, g[1])
        push!(y, g[2])
        push!(psi, abs(scfres.ψ[kidx][gi, bandidx]))
    end
end
Gxmax = maximum(abs.(x))
Gymax = maximum(abs.(y))

z = zeros(Gxmax*2+1, Gymax*2+1)
for (gi, g) in enumerate(scfres.basis.kpoints[kidx].G_vectors)
    if g[3]==Gz
        z[g[1]+Gxmax+1, g[2]+Gymax+1] =abs(scfres.ψ[kidx][gi, bandidx])
    end
end

data = PlotlyJS.contour(x=[i for i in -Gxmax:Gxmax], y=[i for i in -Gymax:Gymax], z=z)
PlotlyJS.plot(data)
# pyplot()
# Plots.contour(x, y, psi, xlabel="Gx", ylabel="Gy", zlabel="|ψ|", legend=false)
# using PlotlyJS
# trace = PlotlyJS.surface(x=x, y=y, z=psi)
# layout = Layout(title="Mt. Bruno Elevation", autosize=false, width=500,
#                     height=500)
# PlotlyJS.plot(trace, layout)