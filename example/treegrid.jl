#--- test tree divided BZ ---

using Random
using LinearAlgebra
using Plots

struct MeshNod
    depth::Int
    pos::Vector{Int}
    sons::Vector{MeshNod}
end

struct TreeMesh
    dim::Int
    root::MeshNod
    latticevectors::Matrix{Float64}
    gridpoints::Matrix{Float64}
end

function TreeMesh(latticevectors::Matrix, weightcap::Float64; depthcap = 5)
    
end

@inline function density(K::AbstractVector)
    me = 0.5
    T = 0.01

    μ = 1.0
    k = norm(K)
    ϵ = k^2 / (2me)
    dϵdk = k / me

    # μ = 0.5
    # ϵ = -(cos(π*K[1]) + cos(π*K[2])) / (2me)
    # dϵdk = sqrt(sin(π*K[1])^2 + sin(π*K[2])^2) / (2me)

    return (exp((ϵ - μ) / T) / T)/(exp((ϵ - μ) / T) + 1.0)^2 * dϵdk
end

function area(gridpoints)
    a = norm(gridpoints[2,:]-gridpoints[1,:])
    b = norm(gridpoints[3,:]-gridpoints[1,:])
    return a*b
end

area(sn::SquareNod) = area(sn.gridpoints)

function weight(gridpoints; N=100)
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
gpcontainer = [K00, K01, K10, K11]

sn = SquareNod!(gridpoints, weight(gridpoints)/500, gpcontainer; treelevelcap = 10)

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

print(sn)
# print(gpcontainer)
println("length:",length(gpcontainer))

X, Y = zeros(Float64, length(gpcontainer)), zeros(Float64, length(gpcontainer))
for (Ki, K) in enumerate(gpcontainer)
    X[Ki], Y[Ki] = K[1], K[2]
end

p = scatter(X, Y)
display(p)
readline()

