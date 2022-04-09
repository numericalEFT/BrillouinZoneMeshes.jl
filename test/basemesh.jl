@testset "Base Mesh" begin

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

    @testset "Uniform Mesh" begin
        δ = 1e-6

        # 2D
        shift = [[0, δ], [0, -δ], [δ, 0], [-δ, 0]]
        origin = [0.0 0.0]
        latvec = [2 0;1 sqrt(3)]'
        # latvec = [1 0; 0 1]'

        umesh = UniformMesh{2, 3}(origin, latvec)

        for i in 1:size(umesh)
            for j in 1:length(shift)
                # println("i=$i, point:$(umesh[i]), shift:$(shift[j])")
                @test i == floor(umesh, umesh[i] + shift[j])
            end
        end

        # 3D
        shift = [[0, 0, δ], [0, 0, -δ], [δ, 0, 0], [-δ, 0, 0], [0, δ, 0], [0, -δ, 0]]
        origin = [0.0 0.0 0.0]
        latvec = [1.0 0 0; 0 1.0 0; 0 0 1.0]'
        

        umesh = UniformMesh{3, 3}(origin, latvec)

        for i in 1:size(umesh)
            for j in 1:length(shift)
                # println("i=$i, point:$(umesh[i]), shift:$(shift[j])")
                @test i == floor(umesh, umesh[i] + shift[j])
            end
        end

    end

end
