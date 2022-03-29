using TightBinding
using LinearAlgebra
using Plots
la = set_Lattice(2,[[1,0],[0,1]])
add_atoms!(la,[0,0])
add_atoms!(la,[0,0])
add_atoms!(la,[0,0])
add_atoms!(la,[0,0])
add_atoms!(la,[0,0])

tmat = [
-0.7    0 -0.4  0.2 -0.1
-0.8    0    0    0    0
 0.8 -1.5    0    0 -0.3
   0  1.7    0    0 -0.1
-3.0    0    0 -0.2    0
-2.1  1.5    0    0    0
 1.3    0  0.2 -0.2    0
 1.7    0    0  0.2    0
-2.5  1.4    0    0    0
-2.1  3.3    0 -0.3  0.7
 1.7  0.2    0  0.2    0
 2.5    0    0  0.3    0
 1.6  1.2 -0.3 -0.3 -0.3
   0    0    0 -0.1    0
 3.1 -0.7 -0.2    0    0
]
tmat = 0.1.*tmat
imap = zeros(Int64,5,5)
count = 0
for μ=1:5
    for ν=μ:5
        global count += 1
        imap[μ,ν] = count
    end
end
Is = [1,-1,-1,1,1,1,1,-1,-1,1,-1,-1,1,1,1]
σds = [1,-1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1]
tmat_σy = tmat[:,:]
tmat_σy[imap[1,2],:] = -tmat[imap[1,3],:]
tmat_σy[imap[1,3],:] = -tmat[imap[1,2],:]
tmat_σy[imap[1,4],:] = -tmat[imap[1,4],:]
tmat_σy[imap[2,2],:] = tmat[imap[3,3],:]
tmat_σy[imap[2,4],:] = tmat[imap[3,4],:]
tmat_σy[imap[2,5],:] = -tmat[imap[3,5],:]
tmat_σy[imap[3,3],:] = tmat[imap[2,2],:]
tmat_σy[imap[3,4],:] = tmat[imap[2,4],:]
tmat_σy[imap[3,5],:] = -tmat[imap[2,5],:]
tmat_σy[imap[4,5],:] = -tmat[imap[4,5],:]

hoppingmatrix = zeros(Float64,5,5,5,5)
hops = [-2,-1,0,1,2]
hopelements = [[1,0],[1,1],[2,0],[2,1],[2,2]]

for μ = 1:5
    for ν=μ:5
        for ii=1:5
            ihop = hopelements[ii][1]
            jhop = hopelements[ii][2]
            #[a,b],[a,-b],[-a,-b],[-a,b],[b,a],[b,-a],[-b,a],[-b,-a]

            #[a,b]
            i = ihop +3
            j = jhop +3
            hoppingmatrix[μ,ν,i,j]=tmat[imap[μ,ν],ii]
            #[a,-b] = σy*[a,b] [1,1] -> [1,-1]
            if jhop != 0
                i = ihop +3
                j = -jhop +3
                hoppingmatrix[μ,ν,i,j]=tmat_σy[imap[μ,ν],ii]
            end

            if μ != ν
                #[-a,-b] = I*[a,b] [1,1] -> [-1,-1],[1,0]->[-1,0]
                i = -ihop +3
                j = -jhop +3
                hoppingmatrix[μ,ν,i,j]=Is[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                #[-a,b] = I*[a,-b] = I*σy*[a,b]  #[2,0]->[-2,0]
                if jhop != 0
                    i = -ihop +3
                    j = jhop +3
                    hoppingmatrix[μ,ν,i,j]=Is[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                end
            end
            #[b,a],[b,-a],[-b,a],[-b,-a]
            if jhop != ihop
                #[b,a] = σd*[a,b]
                i = jhop +3
                j = ihop +3
                hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                #[-b,a] = σd*σy*[a,b]
                if jhop != 0
                    i = -jhop +3
                    j = ihop +3
                    hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                end

                if μ != ν
                    #[-b,-a] = σd*[-a,-b] = σd*I*[a,b]
                    i = -jhop +3
                    j = -ihop +3
                    hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*Is[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                    #[b,-a] = σd*[-a,b] = σd*I*[a,-b] = σd*I*σy*[a,b]  #[2,0]->[-2,0]
                    if jhop != 0
                        i = jhop +3
                        j = -ihop +3
                        hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*Is[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                    end
                end
            end
        end


    end
end

for μ=1:5
    for ν=μ:5
        for i = 1:5
            ih = hops[i]
            for j = 1:5
                jh = hops[j]
                if hoppingmatrix[μ,ν,i,j] != 0.0                
                    add_hoppings!(la,hoppingmatrix[μ,ν,i,j],μ,ν,[ih,jh])
                end
            end
        end
    end
end

onsite = [10.75,10.96,10.96,11.12,10.62]
set_onsite!(la,onsite)

ham = hamiltonian_k(la) #Construct the Hamiltonian
function dispersion_Fe(k)
    eigen_graphene = eigen(ham(k))
    energy = eigen_graphene.values
    U_matrix = eigen_graphene.vectors
    return energy, U_matrix
end


"""
This function returns the diagonalized green's function G_diag of Fe-based superconductor, together with the U matrix which can reproduce the original Green's function by:
    G = U*G_diag*inv(U)
"""
function green_Fe(k, μ, ω)
    energies, U = dispersion_Fe(k)
    Green = zeros(ComplexF64, length(energies), length(energies))
    for i=1:length(energies)
        Green[i,i]= -1/( im*ω - energies[i] + μ) 
    end
    return Green, U
end

if abspath(PROGRAM_FILE) == @__FILE__
    k = [π/4,π/4]
    energies, U = dispersion_Fe(k)
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
    pls = heatmap(kx, ky, energies[1,:,:])
    #pls = surface!(kx, ky, energies[2,:,:], legend=false, fillalpha=0.5)
    #pls = plot_fermisurface_2D(la, Eshift = -0.2, nk = 1000)
    savefig("Fe5_lowest_band.png")

"""
Test band structure of Fe-based 5 orbital superconductor.
Should be consistent with T. Nomura, J. Phys. Soc. Jpn. 78, 034716 (2009)
"""
    set_μ!(la,10.96) #set the chemical potential
    klines = set_Klines()
    kmin = [0,0]
    kmax = [π,0]
    add_Kpoints!(klines,kmin,kmax,"(0,0)","(pi,0)",nk=nk)

    kmin = [π,0]
    kmax = [π,π]
    add_Kpoints!(klines,kmin,kmax,"(pi,0)","(pi,pi)",nk=nk)

    kmin = [π,π]
    kmax = [0,0]
    add_Kpoints!(klines,kmin,kmax,"(pi,pi)","(0,0)",nk=nk)

    pls = calc_band_plot(klines,la)
    savefig("Fe5band.png")
end
