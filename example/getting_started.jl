using BrillouinZoneMeshes

DIM = 2;
# square lattice
# notice the input matrix is transposed
lattice = Matrix([1.0 0; 0 1]');
# create Brillouin zone
br = BZMeshes.Cell(lattice=lattice);
# uniform mesh
umesh = BZMeshes.UniformBZMesh(cell=br, size=(8, 8));
# symmetry reduce map
mm = MeshMaps.MeshMap(umesh);
# reduced mesh
rmesh = ReducedBZMesh(umesh, mm);

println("length of two meshes:", length(umesh), ",", length(rmesh))

println("find nearest point of [1,1] in rmesh")
idx = AbstractMeshes.locate(rmesh, [1, 1])
println("index: $(idx)")
println("point: $(rmesh[idx])")

println("volume of rmesh:", AbstractMeshes.volume(rmesh), "=4π^2")

data = zeros(Float64, length(rmesh))
for (i, p) in enumerate(rmesh)
    data[i] = cos(p[1] / 2) + cos(p[2] / 2)
end
println("computing integral on [-π,π]×[-π,π]")
println("I=∫(cos(x/2)+cos(y/2))dxdy")
println("I=16π=$(16π)")
println("I_num≈", AbstractMeshes.integrate(data, rmesh))