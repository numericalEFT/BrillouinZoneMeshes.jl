@testset "BaryCheb" begin
    println("Testing BaryCheb")

    @testset "1D BaryCheb" begin
        n = 4
        x, w = BaryCheb.barychebinit(n)
        println(x)
        println(w)

        f(t) = t
        F(t) = 0.5*t^2

        data = f.(x)

        # test interp
        @test BaryCheb.barycheb(n, 0.5, data, w, x) ≈ f(0.5)

        # test integrate
        vmat = BaryCheb.vandermonde(x)
        println("vandermonde:",vmat)
        invmat = inv(transpose(vmat))
        println("invmat:",invmat)
        x1, x2 = -0.4, 0.0
        b = BaryCheb.weightcoef(x2, 1, n) - BaryCheb.weightcoef(x1, 1, n)
        println("b:",b)
        intw = BaryCheb.calcweight(invmat, b)
        println("intw:",intw)
        @test sum(intw .* data) ≈ F(x2) - F(x1)
        @test BaryCheb.chebint(n, x1, x2, data, invmat) ≈ F(x2) - F(x1)
    end

end
