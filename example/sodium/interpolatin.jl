# using Plots
using FFTW


# assume f(x) where x = -N/2:N/2-1
# return f(y) where y = -Ni/2:Ni/2-1
function fourier_interpolate(f::AbstractArray{T,D}, Ni::Int64) where {T,D}
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
    return f_interp
end

if abspath(PROGRAM_FILE) == @__FILE__


    Nx = 30
    Nxfine = 100
    x = collect((-Nx/2:Nx/2-1) / Nx)
    y = z = x
    xfine = collect((-Nxfine/2:Nxfine/2-1) / Nxfine)
    yfine = zfine = xfine

    ##### DIM =1 #########################
    f = exp.(sin.(2π * x))
    fine = exp.(sin.(2π * xfine))

    f_interp, xi = fourier_interpolate(f, Nxfine)
    println("DIM = 1 test: ", maximum(abs.(f_interp - fine)))

    # plot(x, f, label="Original samples", markershape=:circle)
    # plot!(xi, f_interp, label="Interpolated values")

    ##### DIM =2 #########################
    f = zeros(Float64, Nx, Nx)
    fine = zeros(Float64, Nxfine, Nxfine)
    f2(x, y) = exp(sin(2π * x) + cos(2π * y))
    for xi in 1:Nx
        for yi in 1:Nx
            f[xi, yi] = f2(x[xi], y[yi])
        end
    end
    for xi in 1:Nxfine
        for yi in 1:Nxfine
            fine[xi, yi] = f2(xfine[xi], yfine[yi])
        end
    end

    f_interp, xi = fourier_interpolate(f, Nxfine)
    println("DIM = 2 test: ", maximum(abs.(f_interp - fine)))

    # plot(x, f[1, :], label="Original samples", markershape=:circle)
    # plot!(xfine, f_interp[1, :], label="Interpolated values")

    ##### DIM =3 #########################
    f = zeros(Float64, Nx, Nx, Nx)
    fine = zeros(Float64, Nxfine, Nxfine, Nxfine)
    f3(x, y, z) = exp(sin(2π * x) + cos(2π * y) + sin(2π * z))
    for xi in 1:Nx
        for yi in 1:Nx
            for zi in 1:Nx
                f[xi, yi, zi] = f3(x[xi], y[yi], z[zi])
            end
        end
    end
    for xi in 1:Nxfine
        for yi in 1:Nxfine
            for zi in 1:Nxfine
                fine[xi, yi, zi] = f3(xfine[xi], yfine[yi], zfine[zi])
            end
        end
    end

    f_interp, xi = fourier_interpolate(f, Nxfine)
    println("DIM = 3 test: ", maximum(abs.(f_interp - fine)))

    # plot(x, f[1, 1, :], label="Original samples", markershape=:circle)
    # plot!(xfine, f_interp[1, 1, :], label="Interpolated values")

end
