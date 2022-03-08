
#--- test tree divided BZ ---

using Random
using LinearAlgebra
using Plots

struct SquareNod
    gridpoints::Matrix{Float64}
    sons::Vector{SquareNod}
end

function minarea(sn::SquareNod)
    return (isempty(sn.sons)) ? area(sn) : min([minarea(ssn) for ssn in sn.sons]...)
end

function SquareNod!(gridpoints, weightcap::Float64, gpcontainer; treelevelcap = 5)

    K00, K01, K10, K11 = Matrix(gridpoints[1,:]'), Matrix(gridpoints[2,:]'), Matrix(gridpoints[3,:]'), Matrix(gridpoints[4,:]')
    K0h = 0.5 .* (K00 .+ K01)
    Kh0 = 0.5 .* (K00 .+ K10)
    K1h = 0.5 .* (K10 .+ K11)
    Kh1 = 0.5 .* (K01 .+ K11)
    Khh = 0.5 .* (K00 .+ K11)

    push!(gpcontainer, K0h)
    push!(gpcontainer, Kh0)
    # push!(gpcontainer, K1h)
    # push!(gpcontainer, Kh1)
    push!(gpcontainer, Khh)

    sons = []
    if treelevelcap > 1 && weight(gridpoints) > weightcap
        sn00 = SquareNod!(Matrix([K00; K0h; Kh0; Khh]), weightcap, gpcontainer;treelevelcap= treelevelcap-1)
        sn01 = SquareNod!(Matrix([K0h; K01; Khh; Kh1]), weightcap, gpcontainer;treelevelcap= treelevelcap-1)
        sn10 = SquareNod!(Matrix([Kh0; Khh; K10; K1h]), weightcap, gpcontainer;treelevelcap= treelevelcap-1)
        sn11 = SquareNod!(Matrix([Khh; Kh1; K1h; K11]), weightcap, gpcontainer;treelevelcap= treelevelcap-1)

        push!(sons, sn00)
        push!(sons, sn01)
        push!(sons, sn10)
        push!(sons, sn11)
    end

    return SquareNod(gridpoints, sons)
end

function SquareNod(gridpoints, weightcap::Float64; treelevelcap = 5)

    K00, K01, K10, K11 = Matrix(gridpoints[1,:]'), Matrix(gridpoints[2,:]'), Matrix(gridpoints[3,:]'), Matrix(gridpoints[4,:]')
    K0h = 0.5 .* (K00 .+ K01)
    Kh0 = 0.5 .* (K00 .+ K10)
    K1h = 0.5 .* (K10 .+ K11)
    Kh1 = 0.5 .* (K01 .+ K11)
    Khh = 0.5 .* (K00 .+ K11)

    sons = []
    if treelevelcap > 1 && weight(gridpoints) > weightcap
        sn00 = SquareNod(Matrix([K00; K0h; Kh0; Khh]), weightcap;treelevelcap= treelevelcap-1)
        sn01 = SquareNod(Matrix([K0h; K01; Khh; Kh1]), weightcap;treelevelcap= treelevelcap-1)
        sn10 = SquareNod(Matrix([Kh0; Khh; K10; K1h]), weightcap;treelevelcap= treelevelcap-1)
        sn11 = SquareNod(Matrix([Khh; Kh1; K1h; K11]), weightcap;treelevelcap= treelevelcap-1)

        push!(sons, sn00)
        push!(sons, sn01)
        push!(sons, sn10)
        push!(sons, sn11)
    end

    return SquareNod(gridpoints, sons)
end

@inline function density(K::AbstractVector)
    me = 0.5
    T = 0.01

    # μ = 1.0
    # k = norm(K)
    # ϵ = k^2 / (2me)
    # dϵdk = k / me

    μ = 0.5
    ϵ = -(cos(π*K[1]) + cos(π*K[2])) / (2me)
    dϵdk = sqrt(sin(π*K[1])^2 + sin(π*K[2])^2) / (2me)

    return (exp((ϵ - μ) / T) / T)/(exp((ϵ - μ) / T) + 1.0)^2 * dϵdk
end

function area(gridpoints)
    a = norm(gridpoints[2,:]-gridpoints[1,:])
    b = norm(gridpoints[3,:]-gridpoints[1,:])
    return a*b
end

area(sn::SquareNod) = area(sn.gridpoints)

function weight(gridpoints; N=400)
    result = 0.0
    a = area(gridpoints)

    K1, K2 = gridpoints[2,:]-gridpoints[1,:], gridpoints[3,:]-gridpoints[1,:]
    for i in 1:N
        K = gridpoints[1, :] .+ rand() .* K1 .+ rand() .* K2
        result += density(K)
    end

    return result * a
end

weight(sn::SquareNod; N=100) = weight(sn.gridpoints; N=N)

gridpoints = Matrix([-1.0 -1.0 ; -1 1 ; 1 -1 ; 1 1])
K00, K01, K10, K11 = Matrix(gridpoints[1,:]'), Matrix(gridpoints[2,:]'), Matrix(gridpoints[3,:]'), Matrix(gridpoints[4,:]')
gpcontainer = [K00,]

sn = SquareNod!(gridpoints, weight(gridpoints)/128^2, gpcontainer; treelevelcap = 12)

println(area(sn))
println(weight(sn))

function Base.print(sn::SquareNod; treelevel = 1)
    for i in 1:treelevel
        print(" ")
    end
    print("area=$(area(sn)),")
    print("weight=$(weight(sn)),")
    print("gridpoints:",sn.gridpoints)
    print("\n")

    for ssn in sn.sons
        print(ssn, ;treelevel = treelevel + 1)
    end
end

# print(sn)

println("minarea = $(minarea(sn))")
# print(gpcontainer)
println("length:",length(gpcontainer))

println("compress:", minarea(sn)*length(gpcontainer)/4)

X, Y = zeros(Float64, length(gpcontainer)), zeros(Float64, length(gpcontainer))
for (Ki, K) in enumerate(gpcontainer)
    X[Ki], Y[Ki] = K[1], K[2]
end

# p = plot(legend = false, size = (1024, 1024))
# scatter!(p, X, Y, marker = :cross, markersize = 2)
# savefig(p, "run/grid.pdf")
# display(p)
# readline()

