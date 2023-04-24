
using BrillouinZoneMeshes
using Test

@testset "Type Stable?" begin
    DIM = 2
    N1, N2 = 3, 5
    lattice = Matrix([1/N1/2 0; 0 1.0/N2/2]') .* 2π
    # so that bzmesh[i,j] = (2i-1,2j-1)
    cell = BZMeshes.Cell(lattice=lattice)
    mesh = BaseMesh.UMesh(br=cell, origin=ones(DIM) ./ 2, size=(N1, N2), shift=zeros(DIM))

    @inferred mesh[1]
end