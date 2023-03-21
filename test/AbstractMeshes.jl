@testset "AbstractMeshes" begin
    using BrillouinZoneMeshes.AbstractMeshes

    function test_func_not_implemented(func, obj)
        # if a func required is not implemented for obj
        # an error occur
        try
            func(obj)
        catch e
            @test e isa ErrorException
        end
    end

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
    test_func_not_implemented(println, notamesh)

    test_func_not_implemented(x -> getindex(x, 1), notamesh)
    test_func_not_implemented(x -> getindex(x, 1, 2, 3), notamesh)
    test_func_not_implemented(x -> getindex(x, FracCoords, 1), notamesh)
    test_func_not_implemented(x -> getindex(x, FracCoords, 1, 2, 3), notamesh)

    test_func_not_implemented(x -> locate(x, 1), notamesh)
    test_func_not_implemented(x -> volume(x, 1), notamesh)
    test_func_not_implemented(volume, notamesh)

    test_func_not_implemented(lattice_vector, notamesh)
    test_func_not_implemented(inv_lattice_vector, notamesh)
    test_func_not_implemented(cell_volume, notamesh)

    test_func_not_implemented(x -> integrate([1,], x), notamesh)
    test_func_not_implemented(x -> interp([1,], x, 1), notamesh)

    test_func_not_implemented(x -> interval(x, 1), notamesh)

end