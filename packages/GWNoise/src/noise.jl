"""
    GaussianNoise(curve, domain::TimeDomain;
                  minimum_frequency=first(curve.frequencies),
                  maximum_frequency=min(last(curve.frequencies), sampling_frequency(domain)/2))

Prepare a finite stationary Gaussian noise segment. The domain sample count and
sampling frequency (Hz) determine the duration and one-sided FFT grid.
The active band must be covered by the curve; other bins have zero power.
DC and, for even sample counts, Nyquist are zero. No taper is applied.
Draws are independent periodic FFT segments, not a continuous noise stream.
"""
struct GaussianNoise{TD <: TimeDomain, FD <: FrequencyDomain}
    time_domain::TD
    frequency_domain::FD
    spectrum::Vector{Float64}
end

function GaussianNoise(
        curve::TabulatedPSD, domain::TimeDomain;
        minimum_frequency = first(curve.frequencies),
        maximum_frequency = min(last(curve.frequencies), sampling_frequency(domain) / 2),
    )
    length(channels(domain)) == 1 || throw(ArgumentError("GaussianNoise requires a single-channel domain"))
    fs = sampling_frequency(domain)
    n = nsamples(domain)
    isfinite(fs) && fs > 0 && n >= 2 ||
        throw(ArgumentError("need a positive sampling frequency and at least two samples"))
    first(curve.frequencies) <= minimum_frequency <= maximum_frequency <=
        min(last(curve.frequencies), fs / 2) ||
        throw(ArgumentError("requested band must lie within PSD coverage and Nyquist"))
    fs64 = Float64(fs)
    isfinite(fs64) && fs64 > 0 || throw(ArgumentError("noise sampling frequency must be representable as positive finite Float64"))
    df = fs64 / n
    f = range(0.0; step = df, length = n ÷ 2 + 1)
    s = [minimum_frequency <= x <= maximum_frequency ? psd(x, curve) : 0.0 for x in f]
    s[1] = 0.0
    iseven(n) && (s[end] = 0.0)
    fd = FrequencyDomain(; frequencies = f, epoch = epoch(domain), channels = channels(domain))
    return GaussianNoise(domain, fd, s)
end

"""`frequencies(model)` returns the one-sided FFT frequency range in Hz."""
frequencies(model::GaussianNoise) = frequencies(model.frequency_domain)

"""`psd(model)` returns a copy of the effective PSD, including band and endpoint zeros."""
psd(model::GaussianNoise) = copy(model.spectrum)

"""
    noise_fd([rng], model::GaussianNoise)

Draw complex Fourier strain with the normalization `rfft(strain)/sampling_frequency`.
For each active interior bin, the real and imaginary parts independently have
variance `PSD/(4Δf)`. Each call draws a new realization.
"""
noise_fd(model::GaussianNoise) = noise_fd(default_rng(), model)
function noise_fd(rng::AbstractRNG, model::GaussianNoise)
    out = Vector{ComplexF64}(undef, length(frequencies(model)))
    noise_fd!(rng, out, model)
    return out
end

"""`noise_fd!([rng], out, model)` fills a complex vector with a new FD draw."""
noise_fd!(out::AbstractVector{<:Complex}, model::GaussianNoise) =
    noise_fd!(default_rng(), out, model)

function noise_fd!(rng::AbstractRNG, out::AbstractVector{<:Complex}, model::GaussianNoise)
    df = frequency_spacing(model.frequency_domain)
    for i in eachindex(out, model.spectrum)
        scale = sqrt(model.spectrum[i] / df) / 2
        out[i] = scale * complex(randn(rng), randn(rng))
    end
    return out
end

"""
    noise_td([rng], model::GaussianNoise)

Draw real time-domain strain using `to_timedomain(model.time_domain, noise_fd(rng, model))`.
To obtain matching FD and TD data, draw once with `noise_fd` and transform it.
"""
noise_td(model::GaussianNoise) = noise_td(default_rng(), model)
noise_td(rng::AbstractRNG, model::GaussianNoise) =
    to_timedomain(model.time_domain, noise_fd(rng, model))

"""Validate a stationary noise model against sampling and channels, ignoring epoch."""
function check_compatible(d::TimeDomain, model::GaussianNoise)
    nsamples(d) == nsamples(model.time_domain) || throw(DimensionMismatch("noise sample count differs"))
    fs = Float64(sampling_frequency(model.time_domain))
    check_coordinate(sampling_frequency(d), fs, fs, "sampling frequencies")
    channels(d) == channels(model.time_domain) || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
function check_compatible(d::FrequencyDomain, model::GaussianNoise)
    check_same_frequency_grid(frequencies(d), frequencies(model.frequency_domain))
    channels(d) == channels(model.frequency_domain) || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
