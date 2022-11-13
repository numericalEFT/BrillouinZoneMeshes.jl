using Plots
using FFTW

Nx = 30;
x = collect((-Nx/2:Nx/2-1) / Nx);

f = exp.(sin.(2π * x));
# ft = fftshift(fft(f))

# Npad = floor(Int64, Ni / 2 - N / 2)
# ft_pad = [zeros(Npad); ft; zeros(Npad)];
# f_interp = real(ifft(fftshift(ft_pad))) * Ni / N;


# assume f(x) where x = -N/2:N/2-1
# return f(y) where y = -Ni/2:Ni/2-1
function interpolate(f::AbstractVector, Ni::Int64)
    N = length(f)
    Npad = floor(Int64, Ni / 2 - N / 2)
    ft = fftshift(fft(fftshift(f))) # first shift t0 0: N-1, then fft
    ft_pad = [zeros(Npad); ft; zeros(Npad)]
    f_interp = real(fftshift(ifft(fftshift(ft_pad)))) * Ni / N
    xi = collect((-Ni/2:Ni/2-1) / Ni)
    return f_interp, xi
end

function interpolate(f::AbstractArray, Ni::Int64, dim::Int)
    N = size(f, dim)
    Npad = floor(Int64, Ni / 2 - N / 2)
    sf = fftshift(f, dim) # first shift t0 0: N-1
    fsf = fft!(sf, dim) # then fft
    ft = fftshift!(fsf, dim) # shift to N/2: N/2-1
    # TODO: doesn't work for high dimension for now
    # TODO: how to write a padding for generic dimension?
    ft_pad = [zeros(Npad); ft; zeros(Npad)] # pad to Ni/2: Ni/2-1
    sft_pad = fftshift!(ft_pad, dim) # shift to 0: Ni-1
    isft_pad = ifft!(sft_pad, dim) # ifft
    f_interp = real(fftshift(isft_pad)) * Ni / N #shift to Ni/2: Ni/2-1
    return f_interp
end

f_interp, xi = interpolate(f, 300)
plot(x, f, label="Original samples", markershape=:circle)
plot!(xi, f_interp, label="Interpolated values")