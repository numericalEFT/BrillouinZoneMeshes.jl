@testset "CompositeMeshes.jl" begin
    using BrillouinZoneMeshes.CompositeGrids
    using BrillouinZoneMeshes.BaseMesh
    using BrillouinZoneMeshes.AbstractMeshes
    using BrillouinZoneMeshes.CompositeMeshes

    rng = MersenneTwister(1234)

    x = SimpleGrid.Uniform([0.0, π], 3)
    y = SimpleGrid.Uniform([0.0, π], 3)

    dpm = DirectProdMesh(x, y)
    cm = CompositeMesh(dpm, 3)

    vol = 0.0
    for (pi, p) in enumerate(cm)
        i, j = AbstractMeshes._ind2inds(size(cm), pi)
        @test pi == AbstractMeshes.locate(cm, p)
        vol += AbstractMeshes.volume(cm, pi)
        #println(p)
    end
    @test vol ≈ AbstractMeshes.volume(cm)

    f(x) = sin(x[1]) + cos(x[2])

    data = zeros(Float64, size(cm))

    for i in 1:length(cm)
        data[i] = f(cm[i])
    end

    ## interpolate
    testN = 3
    xlist = rand(rng, testN) * π
    ylist = rand(rng, testN) * π
    for x in xlist
        for y in ylist
            @test isapprox(f([x, y]), interp(data, cm, [x, y]), rtol=4e-2)
        end
    end

    ## integrate
    integral = integrate(data, cm)
    # println("integral=$(integral)")
    @test isapprox(integral, 2π, rtol=1e-3)

end