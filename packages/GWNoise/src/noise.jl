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
struct GaussianNoise{F <: AbstractRange, C <: Tuple}
    sampling_frequency::Float64
    nsamples::Int
    frequencies::F
    spectrum::Vector{Float64}
    scales::Vector{Float64}
    channels::C
end

function GaussianNoise(
        curve::TabulatedPSD, domain::TimeDomain;
        minimum_frequency = first(curve.frequencies),
        maximum_frequency = min(last(curve.frequencies), GWTimeFrequencyDomains.sampling_frequency(domain) / 2),
    )
    length(channels(domain)) == 1 || throw(ArgumentError("GaussianNoise requires a single-channel domain"))
    sampling_frequency = GWTimeFrequencyDomains.sampling_frequency(domain)
    nsamples = GWTimeFrequencyDomains.nsamples(domain)
    isfinite(sampling_frequency) && sampling_frequency > 0 && nsamples >= 2 ||
        throw(ArgumentError("need a positive sampling frequency and at least two samples"))
    first(curve.frequencies) <= minimum_frequency <= maximum_frequency <=
        min(last(curve.frequencies), sampling_frequency / 2) ||
        throw(ArgumentError("requested band must lie within PSD coverage and Nyquist"))
    fs = Float64(sampling_frequency)
    isfinite(fs) && fs > 0 || throw(ArgumentError("noise sampling frequency must be representable as positive finite Float64"))
    df = fs / nsamples
    f = range(0.0; step = df, length = nsamples ÷ 2 + 1)
    s = [minimum_frequency <= x <= maximum_frequency ? psd(x, curve) : 0.0 for x in f]
    s[1] = 0.0
    iseven(nsamples) && (s[end] = 0.0)
    scales = sqrt.(s ./ df) ./ 2
    return GaussianNoise(fs, Int(nsamples), f, s, scales, channels(domain))
end

"""`frequencies(model)` returns the one-sided FFT frequency range in Hz."""
frequencies(model::GaussianNoise) = model.frequencies

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
    out = Vector{ComplexF64}(undef, length(model.frequencies))
    noise_fd!(rng, out, model)
    return out
end

"""`noise_fd!([rng], out, model)` fills a complex vector with a new FD draw."""
noise_fd!(out::AbstractVector{<:Complex}, model::GaussianNoise) =
    noise_fd!(default_rng(), out, model)

function noise_fd!(rng::AbstractRNG, out::AbstractVector{<:Complex}, model::GaussianNoise)
    for i in eachindex(out, model.scales)
        out[i] = model.scales[i] * complex(randn(rng), randn(rng))
    end
    return out
end

"""
    noise_td([rng], model::GaussianNoise)

Draw real time-domain strain: `irfft(noise_fd(rng, model), nsamples) * sampling_frequency`.
To obtain matching FD and TD data, draw once with `noise_fd` and transform it.
"""
noise_td(model::GaussianNoise) = noise_td(default_rng(), model)
noise_td(rng::AbstractRNG, model::GaussianNoise) =
    irfft(noise_fd(rng, model), model.nsamples) .* model.sampling_frequency

"""
    noise_td!([rng], out, model)

Fill a real vector with a new TD draw. This convenience method allocates FFT
temporaries; it does not cache scratch buffers or FFT plans in the model.
"""
noise_td!(out::AbstractVector{<:Real}, model::GaussianNoise) =
    noise_td!(default_rng(), out, model)

function noise_td!(rng::AbstractRNG, out::AbstractVector{<:Real}, model::GaussianNoise)
    return copyto!(out, noise_td(rng, model))
end

"""Validate a stationary noise model against sampling and channels, ignoring epoch."""
function check_compatible(d::TimeDomain, model::GaussianNoise)
    nsamples(d) == model.nsamples || throw(DimensionMismatch("noise sample count differs"))
    GWTimeFrequencyDomains.check_coordinate(sampling_frequency(d), model.sampling_frequency, model.sampling_frequency, "sampling frequencies")
    channels(d) == model.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
function check_compatible(d::FrequencyDomain, model::GaussianNoise)
    check_same_frequency_grid(frequencies(d), frequencies(model))
    channels(d) == model.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
