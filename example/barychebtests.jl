using CodeTracking

using BrillouinZoneMeshes
using BrillouinZoneMeshes.AbstractTrees, BrillouinZoneMeshes.GridTree, BrillouinZoneMeshes.BaseMesh, BrillouinZoneMeshes.BaryCheb

macro expr2fn(fname, expr, args...)
    fn = quote
        function $(esc(fname))()
            $(esc(expr.args[1]))
        end
    end
    for arg in args
        push!(fn.args[2].args[1].args, esc(arg))
    end
    return fn
end

function test1D(f, F; Ni=3, Nf=12)

    for n in Ni:Nf
        bc = BaryCheb.BaryCheb1D(n)

        # data1 = f.(bc.x)
        data1 = [f(xi) for xi in bc.x]
        analytic1 = F(1) - F(-1)
        numeric1 = BaryCheb.integrate1D(data1, bc)
        println("n=$n, $analytic1 <-> $numeric1, diff:$(abs(analytic1-numeric1))")
        # println("n=$n")
        # println("f(x)=cos(x), $analytic1 <-> $numeric1, diff:$(abs(analytic1-numeric1)); f(x)=ln(x+x0), $analytic2 <-> $numeric2, diff:$(abs(analytic2-numeric2))")
        # println("f(x)=ln(x+1), $analytic2 <-> $numeric2, diff:$(abs(analytic2-numeric2))")
    end
end

# test1D(:(cos(x)), :(sin(x)))
println("cos <-> sin")
test1D(x -> (cos(x)), x -> (sin(x)))

x0 = 1.5
println("ln(x+x0) <-> (x+x0)ln(x+x0)-x")
test1D(x -> (log(x + x0)), x -> ((x + x0) * log(x + x0) - x))

println("x^4 <-> x^5/5")
test1D(x -> x^4, x -> x^5 / 5.0)
