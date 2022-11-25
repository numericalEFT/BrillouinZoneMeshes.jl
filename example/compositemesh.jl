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
T = 0.001

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

N = 8
bound = [-π, π]
theta = SimpleGrid.Uniform(bound, N; isperiodic=true)
# bzmesh = BZMeshes.PolarMesh(dispersion=dispersion, anglemesh=theta, cell=br,
#     kmax=π, Nloggrid=5, Nbasegrid=4, minterval=0.001)
pm = BZMeshes.CompositePolarMesh(dispersion=dispersion,
    anglemesh=theta, cell=br, basegridtype=:uniform,
    kmax=π * sqrt(2.4), Nloggrid=8, Nbasegrid=3, minterval=0.1T, N=4)

pm2 = BZMeshes.CompositePolarMesh(dispersion=dispersion,
    anglemesh=theta, cell=br, basegridtype=:uniform,
    kmax=π * sqrt(2.4), Nloggrid=10, Nbasegrid=6, minterval=0.01T, N=6)


println("N=$(length(pm)), ", int_mesh(pm))
println("N=$(length(pm2)), ", int_mesh(pm2))

Nuni = 1000
um = BZMeshes.UniformBZMesh(cell=br, size=(Nuni, Nuni))
um2 = BZMeshes.UniformBZMesh(cell=br, size=(2Nuni, 2Nuni))
println("N=$(length(um)), ", int_mesh(um))
println("N=$(length(um2)), ", int_mesh(um2))