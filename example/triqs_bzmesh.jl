using PythonCall
using SpaceGrid

gf = pyimport("triqs.gf")
lat = pyimport("triqs.lattice")
tb = pyimport("triqs.lattice.tight_binding")
np = pyimport("numpy")
tpl = pyimport("triqs.plot.mpl_interface")
plt = pyimport("matplotlib.pyplot")

BL = lat.BravaisLattice(units=((1, 0, 0), (0, 1, 0))) #square lattice
BZ = lat.BrillouinZone(BL)
nk = 5
mk = gf.MeshBrillouinZone(BZ, nk)

println(BZ.units)
latvec = pyconvert(Array, BZ.units)[1:2, 1:2]
println(latvec)

umesh = SpaceGrid.BaseMesh.UniformMesh{2, nk, SpaceGrid.BaseMesh.EdgedMesh}([0.0,0.0], latvec)

for (ip, p) in enumerate(mk)
    println(p, umesh[ip])
end
