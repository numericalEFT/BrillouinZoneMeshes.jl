w = zeros(4)
et = [1.0, 2.0, 3.0, 4.0]
ef = 1.5

ccall((:intweight1t_, "./tetra.so"), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}), w, et, ef)

println("weight: $w")