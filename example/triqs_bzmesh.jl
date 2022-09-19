using PythonCall
using BZMeshes

gf = pyimport("triqs.gf")
lat = pyimport("triqs.lattice")
tb = pyimport("triqs.lattice.tight_binding")
np = pyimport("numpy")
tpl = pyimport("triqs.plot.mpl_interface")
plt = pyimport("matplotlib.pyplot")

BL = lat.BravaisLattice(units=((1, 0, 0), (0, 1, 0))) #square lattice
BZ = lat.BrillouinZone(BL)
nk = 8
mk = gf.MeshBrillouinZone(BZ, nk)
miw = gf.MeshImFreq(beta=1.0, S="Fermion", n_max=100) #grid number : 201
mprod = gf.MeshProduct(mk, miw)

println(BZ.units)
latvec = pyconvert(Array, BZ.units)[1:2, 1:2]
println(latvec)

umesh = BZMeshes.BaseMesh.UniformMesh{2,nk,BZMeshes.BaseMesh.EdgedMesh}([0.0, 0.0], latvec)

for (ip, p) in enumerate(mk)
    println(p, umesh[ip])
end

G_w = gf.GfImFreq(mesh=miw, target_shape=[1, 1]) #G_w.data.shape will be [201, 1, 1]
G_k_w = gf.GfImFreq(mesh=mprod, target_shape=[1, 1]) #G_k_w.data.shape will be [400, 201, 1, 1

t = 1.0
U = 4.0

G_k_w.data.fill(0.0)
for (ik, k) in enumerate(G_k_w.mesh[0])
    G_w << gf.inverse(gf.iOmega_n - 2 * t * (np.cos(k[0]) + np.cos(k[1])))
    # G_k_w.data[ik-1, pyslice(nothing), 0, 0] = G_w.data[pyslice(nothing), 0, 0] #pyslice(nothing) == :, maybe there is a better way to do this
    G_k_w.data[ik-1, pyslice(0, nk^2), 0, 0] = G_w.data[pyslice(0, nk^2), 0, 0] #pyslice(nothing) == :, maybe there is a better way to do this
end

data0 = pyconvert(Array, G_k_w.data[pyslice(0, nk^2), 0, 0, 0])

data1 = zeros(ComplexF64, (nk, nk))

for p in mk
    inds = pyconvert(Array, p.index)[1:2] .+ 1
    pval = pyconvert(Array, p.value)
    println(pval, umesh[inds...])
    data1[inds...] = pyconvert(ComplexF64, G_k_w.data[p.linear_index, 0, 0, 0])
end

println(data0[1:8])
println(data1[1:8])

for i in 1:nk^2
    println(data0[i] ≈ data1[i])
end
