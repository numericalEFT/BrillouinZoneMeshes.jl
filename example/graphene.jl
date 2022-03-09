using TightBinding
using Plots
using LinearAlgebra
#Primitive vectors
a1 = [sqrt(3)/2,1/2]
a2= [0,1]
#set lattice
la = set_Lattice(2,[a1,a2])
#add atoms
add_atoms!(la,[1/3,1/3])
add_atoms!(la,[2/3,2/3])
show_neighbors(la)
t = 1.0
add_hoppings!(la,-t,1,2,[1/3,1/3])
add_hoppings!(la,-t,1,2,[-2/3,1/3])
add_hoppings!(la,-t,1,2,[1/3,-2/3])
plot_lattice_2d(la)
nk = 100 #numer ob meshes. nk^d meshes are used. d is a dimension.
plot_DOS(la, nk)
#show the band structure
klines = set_Klines()
kmin = [0,0]
kmax = [2π/sqrt(3),0]
add_Kpoints!(klines,kmin,kmax,"G","K")

kmin = [2π/sqrt(3),0]
kmax = [2π/sqrt(3),2π/3]
add_Kpoints!(klines,kmin,kmax,"K","M")

kmin = [2π/sqrt(3),2π/3]
kmax = [0,0]
add_Kpoints!(klines,kmin,kmax,"M","G")
calc_band_plot(klines,la)

ham = hamiltonian_k(la) #Construct the Hamiltonian


function dispersion_graphene(k)
    eigen_graphene = eigen(ham(k))
    energy = eigen_graphene.values
    U_matrix = eigen_graphene.vectors
    return energy, U_matrix
end


"""
This function returns the diagonalized green's function G_diag of graphene, together with the U matrix which can reproduce the original Green's function by:
    G = U*G_diag*inv(U)
"""
function green_graphene(k, μ, ω)
    energies, U = dispersion_graphene(k)
    Green = zeros(ComplexF64, length(energies), length(energies))
    for i=1:length(energies)
        Green[i,i]= -1/( im*ω - energies[i] + μ) 
    end
    return Green, U
end

if abspath(PROGRAM_FILE) == @__FILE__
    k = [π/4,π/4]
    energies, U = dispersion_graphene(k)
    #U = transpose(U)
    Ematrix = zeros(Float64, length(energies), length(energies))
    for i=1:length(energies)
        Ematrix[i,i]=energies[i]
    end
    dim = la.dim
    n = la.numatoms
    print("Error of matrix diagonalization is ", maximum(abs.(U*Ematrix*inv(U)-ham(k))),"\n")
    nk = 100
    energies = zeros(Float64,n,nk,nk)
    kx = range(-2π, stop = 2π, length = nk)
    ky = range(-2π, stop = 2π, length = nk)
    for i1=1:nk
        for i2=1:nk
            energy = dispersion(dim,n,ham,[kx[i1],ky[i2]])
            for j=1:n
                energies[j,i1,i2] = energy[j]
            end

        end
    end
    #pls = heatmap(kx, ky, energies[1,:,:])
    pls = surface(kx, ky, energies[1,:,:], legend=false, fillalpha=0.5)
    pls = surface!(kx, ky, energies[2,:,:], legend=false, fillalpha=0.5)
    #pls = plot_fermisurface_2D(la, Eshift = -0.2, nk = 1000)
    savefig("graphene_band.png")
end
