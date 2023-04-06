@testset "AbstractMeshes" begin
    using BrillouinZoneMeshes.AbstractMeshes

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

    # indexing
    @test mesh[end] == mesh[length(mesh)]
    @test mesh[begin] == mesh[1]

    # tools
    @test AbstractMeshes._inds2ind(size(mesh), (2, 2)) == 5
    @test AbstractMeshes._ind2inds(size(mesh), 5) == [2, 2]

    # test error thrown from funcs not implemented
    struct NotAMesh{T,DIM} <: AbstractMesh{T,DIM} end
    notamesh = NotAMesh{Float64,3}()
    @test_throws ErrorException println(notamesh)

    @test_throws ErrorException notamesh[1]
    @test_throws ErrorException notamesh[1, 2, 3]
    @test_throws ErrorException notamesh[FracCoords, 1]
    @test_throws ErrorException notamesh[FracCoords, 1, 2, 3]

    @test_throws ErrorException locate(notamesh, 1)
    @test_throws ErrorException volume(notamesh, 1)
    @test_throws ErrorException volume(notamesh)

    @test_throws ErrorException lattice_vector(notamesh)
    @test_throws ErrorException inv_lattice_vector(notamesh)
    @test_throws ErrorException cell_volume(notamesh)

    @test_throws ErrorException integrate([1,], notamesh)
    @test_throws ErrorException interp([1,], notamesh, 1)

    @test_throws ErrorException interval(notamesh, 1)

end