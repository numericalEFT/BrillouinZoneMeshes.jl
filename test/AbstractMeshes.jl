@testset "AbstractMeshes" begin

    # create a random concrete mesh
    DIM = 2
    N1, N2 = 3, 5
    lattice = Matrix([1/N1/2 0; 0 1.0/N2/2]') .* 2π
    cell = BZMeshes.Cell(lattice=lattice)
    mesh = BaseMesh.UMesh(br=cell, origin=ones(DIM) ./ 2, size=(N1, N2), shift=zeros(DIM))

    # type
    @test eltype(mesh) == BrillouinZoneMeshes.SVector{Float64,DIM}

    # size
    @test length(mesh) == N1 * N2
    @test size(mesh) == (N1, N2)
    @test size(mesh, 1) == N1
    @test ndims(mesh) == DIM

    # tools
    @test AbstractMeshes._inds2ind(size(mesh), (2, 2)) == 5
    @test AbstractMeshes._ind2inds(size(mesh), 5) == [2, 2]
end