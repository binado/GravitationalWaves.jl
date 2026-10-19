# Local adapters: do not add methods to FFTW or ForwardDiff.
"""FFTW transform of real values, preserving first-order ForwardDiff partials."""
real_fft(x::AbstractVector{<:AbstractFloat}) = rfft(x)

function real_fft(x::AbstractVector{ForwardDiff.Dual{Tag, V, P}}) where {Tag, V <: AbstractFloat, P}
    primal = rfft(ForwardDiff.value.(x))
    partials = ntuple(j -> rfft(map(v -> ForwardDiff.partials(v)[j], x)), P)
    return map(eachindex(primal)) do k
        re = ForwardDiff.Dual{Tag}(real(primal[k]), ntuple(j -> real(partials[j][k]), P))
        im = ForwardDiff.Dual{Tag}(imag(primal[k]), ntuple(j -> imag(partials[j][k]), P))
        complex(re, im)
    end
end

"""Physical Δt*rfft, corrected for the nonzero absolute sample start time."""
function physical_spectrum(x, dt, t_start)
    N = length(x)
    return real_fft(x) .* [dt * cis(-2π * k * t_start / (N * dt)) for k in 0:div(N, 2)]
end

quintic_step(x) = x <= 0 ? zero(x) : x >= 1 ? one(x) : x^3 * (10 + x * (-15 + 6x))
