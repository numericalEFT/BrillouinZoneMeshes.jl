using Plots
using FFTW

dim = 2

Nx = 30;
Nxfine = 60;
x = collect((-Nx/2:Nx/2-1) / Nx);
y = z = x
xfine = collect((-Nxfine/2:Nxfine/2-1) / Nxfine);
yfine = zfine = xfine

if dim == 1
    f = exp.(sin.(2π * x))
    fine = exp.(sin.(2π * xfine))

elseif dim == 2
    f = zeros(Float64, Nx, Nx)
    fine = zeros(Float64, Nxfine, Nxfine)
    for xi in 1:Nx
        for yi in 1:Nx
            f[xi, yi] = exp(sin(2π * x[xi]) + cos(2π * y[yi]))
        end
    end
    for xi in 1:Nxfine
        for yi in 1:Nxfine
            fine[xi, yi] = exp(sin(2π * xfine[xi]) + cos(2π * yfine[yi]))
        end
    end
else
    error("not implemented!")
end

# ft = fftshift(fft(f))

# Npad = floor(Int64, Ni / 2 - N / 2)
# ft_pad = [zeros(Npad); ft; zeros(Npad)];
# f_interp = real(ifft(fftshift(ft_pad))) * Ni / N;


# function interpolate(f::AbstractVector, Ni::Int64)
#     N = length(f)
#     Npad = floor(Int64, Ni / 2 - N / 2)
#     ft = fftshift(fft(fftshift(f))) # first shift t0 0: N-1, then fft
#     ft_pad = [zeros(Npad); ft; zeros(Npad)]
#     f_interp = real(fftshift(ifft(fftshift(ft_pad)))) * Ni / N
#     xi = collect((-Ni/2:Ni/2-1) / Ni)
#     return f_interp, xi
# end

# assume f(x) where x = -N/2:N/2-1
# return f(y) where y = -Ni/2:Ni/2-1
function interpolate(f::AbstractArray{T,D}, Ni::Int64) where {T,D}
    N = size(f, 1)
    @assert all(x -> x == N, size(f))
    dim = [i for i in 1:D]
    println(dim)
    Npad = floor(Int64, Ni / 2 - N / 2)
    sf = fftshift(f, dim) # first shift t0 0: N-1
    fsf = fft(sf, dim) # then fft
    ft = fftshift(fsf, dim) # shift to N/2: N/2-1
    if D == 1
        ft_pad = zeros(eltype(ft), Ni)
        ft_pad[Npad+1:Npad+N] .= ft
        # ft_pad = [zeros(Npad); ft; zeros(Npad)] # pad to Ni/2: Ni/2-1
    elseif D == 2
        ft_pad = zeros(eltype(ft), Ni, Ni)
        ft_pad[Npad+1:Npad+N, Npad+1:Npad+N] .= ft
    elseif D == 3
        ft_pad = zeros(eltype(ft), Ni, Ni, Ni)
        ft_pad[Npad+1:Npad+N, Npad+1:Npad+N, Npad+1:Npad+N] .= ft
    else
        error("not implemented")
    end
    sft_pad = fftshift(ft_pad, dim) # shift to 0: Ni-1
    isft_pad = ifft(sft_pad, dim) # ifft
    f_interp = real(fftshift(isft_pad)) * (Ni / N)^D #shift to Ni/2: Ni/2-1
    xi = collect((-Ni/2:Ni/2-1) / Ni)
    return f_interp, xi
end

if dim == 1
    f_interp, xi = interpolate(f, 300)
    plot(x, f, label="Original samples", markershape=:circle)
    plot!(xi, f_interp, label="Interpolated values")
elseif dim == 2
    # f_interp, xi = interpolate(f, 300)
    # # p = plot(st= :surface, camera = (30, 60), size=(600, 600))
    # # plot!(p, x, y, f, st= :surface)
    # # plot!(p, xi, xi, f_interp, st= :surface)

    f_interp, xi = interpolate(f, Nxfine)
    println(maximum(abs.(f_interp - fine)))

    plot(x, f[1, :], label="Original samples", markershape=:circle)
    plot!(xi, f_interp[1, :], label="Interpolated values")
else
    error("not implemented!")
end