@testset "Tree Grid" begin
    function dispersion(k)
        me = 0.5
        return dot(k, k) / 2me
    end

	  function density(k)
        T = 0.01
        μ = 1.0

        ϵ = dispersion(k) - μ

        # return 1 / (exp((ϵ) / T) + 1.0)
        return (π*T)^2/((π*T)^2 + ϵ^2)
        # return (exp((ϵ) / T) / T)/(exp((ϵ) / T) + 1.0)^2
    end

    latvec = [2 0; 1 sqrt(3)]' .* (2π)
    tg = treegridfromdensity(k->density(k), latvec; atol = 1/2^12, maxdepth = 6, mindepth = 1, N = 2)

    println("size:$(size(tg)),\t length:$(length(tg)),\t efficiency:$(efficiency(tg))")

    @testset "SymMap" begin
        atol = 1e-6
        smap = SymMap(tg, k->density(k); atol = atol)
        println(smap.map)
        println("compress:$(smap.reduced_length/length(smap.map))")

        vals = smap._vals

        # test map
        for i in 1:size(tg)
            @test isapprox(density(tg[i]), vals[smap.map[i]], atol = 2atol)
        end

        # test inv_map
        for i in 1:length(smap.inv_map)
            val = density(tg[smap.inv_map[i][1]])
            for j in 1:length(smap.inv_map[i])
                @test isapprox(val, density(tg[smap.inv_map[i][j]]), atol = 2atol)
            end
        end
    end
end
