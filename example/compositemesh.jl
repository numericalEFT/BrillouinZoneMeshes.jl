using BrillouinZoneMeshes
using BrillouinZoneMeshes.BZMeshes
using BrillouinZoneMeshes.BZMeshes.Coordinates
using BrillouinZoneMeshes.CompositeGrids
using BrillouinZoneMeshes.BaseMesh
using BrillouinZoneMeshes.AbstractMeshes
using BrillouinZoneMeshes.LinearAlgebra
using BrillouinZoneMeshes.StaticArrays

lattice = Matrix([1 0; 0 1.0]')
atoms = [1]
positions = [zeros(2)]
br = BZMeshes.Cell(lattice=lattice, atoms=atoms, positions=positions)

_dispersion(k) = -sum(cos.(k)) + 0.1
# _dispersion(k) = dot(k, k) - π^2 / 4
T = 0.1

function dispersion(k)
    # constrain to 1st bz
    _k = [(abs(p) < π) ? p : π for p in k]
    return _dispersion(_k)
end

function integrand(k)
    #if dot(k, k) <= π^2
    if abs(k[1]) <= π && abs(k[2]) <= π
        return 1 / (dispersion(k)^2 + T^2)
    else
        return 0.0
    end
end

function int_mesh(mesh)
    data = zeros(Float64, size(mesh))
    for i in 1:length(mesh)
        data[i] = integrand(mesh[i])
    end
    return integrate(data, mesh)
end

N = 16
bound = [-π, π]
theta = SimpleGrid.Uniform(bound, N; isperiodic=true)
# bzmesh = BZMeshes.PolarMesh(dispersion=dispersion, anglemesh=theta, cell=br,
#     kmax=π, Nloggrid=5, Nbasegrid=4, minterval=0.001)
pm = BZMeshes.CompositePolarMesh(dispersion=dispersion,
    anglemesh=theta, cell=br, basegridtype=:uniform,
    kmax=π * sqrt(2.1), Nloggrid=12, Nbasegrid=2, minterval=0.1T, N=4)

pm2 = BZMeshes.CompositePolarMesh(dispersion=dispersion,
    anglemesh=theta, cell=br, basegridtype=:uniform,
    kmax=π * sqrt(2.1), Nloggrid=8, Nbasegrid=4, minterval=0.5T, N=6)


println("N=$(length(pm)), ", int_mesh(pm))
println("N=$(length(pm2)), ", int_mesh(pm2))

um = BZMeshes.UniformBZMesh(cell=br, size=(1000, 1000))
um2 = BZMeshes.UniformBZMesh(cell=br, size=(2000, 2000))
println("N=$(length(um)), ", int_mesh(um))
println("N=$(length(um2)), ", int_mesh(um2))