"""Prepared finite stationary Gaussian noise segment, with metadata independent of storage."""
struct GaussianNoise{TM <: NamedTuple, FM <: NamedTuple}
    time_metadata::TM
    frequency_metadata::FM
    nsamples::Int
    spectrum::Vector{Float64}
end

function GaussianNoise(curve::TabulatedPSD, series::TimeSeries; kwargs...)
    return GaussianNoise(curve; nsamples = nsamples(series), metadata(series)..., kwargs...)
end
function GaussianNoise(
        curve::TabulatedPSD; nsamples::Integer, epoch::Real, channels,
        sampling_interval = nothing, sampling_frequency = nothing,
        minimum_frequency = first(curve.frequencies), maximum_frequency = nothing,
    )
    tm = GWSeries.time_metadata(; epoch, channels, sampling_interval, sampling_frequency)
    length(tm.channels) == 1 || throw(ArgumentError("GaussianNoise requires one channel"))
    nsamples >= 2 || throw(ArgumentError("noise needs at least two samples"))
    fs = inv(tm.sampling_interval)
    upper = something(maximum_frequency, min(last(curve.frequencies), fs / 2))
    first(curve.frequencies) <= minimum_frequency <= upper <= min(last(curve.frequencies), fs / 2) ||
        throw(ArgumentError("requested band must lie within PSD coverage and Nyquist"))
    fs64 = Float64(fs)
    isfinite(fs64) && fs64 > 0 || throw(ArgumentError("noise sampling frequency must be representable as positive finite Float64"))
    df = fs64 / nsamples
    f = range(0.0; step = df, length = nsamples ÷ 2 + 1)
    spectrum = [minimum_frequency <= x <= upper ? psd(x, curve) : 0.0 for x in f]
    spectrum[1] = 0.0
    iseven(nsamples) && (spectrum[end] = 0.0)
    fm = (; epoch = tm.epoch, frequency_start = 0.0, frequency_spacing = df, channels = tm.channels, spectrum = :onesided)
    return GaussianNoise(tm, fm, Int(nsamples), spectrum)
end
frequencies(model::GaussianNoise) = range(model.frequency_metadata.frequency_start; step = model.frequency_metadata.frequency_spacing, length = length(model.spectrum))
psd(model::GaussianNoise) = copy(model.spectrum)

"""Draw a one-channel frequency series with interior variance PSD/(4Δf) per component."""
noise_fd(model::GaussianNoise) = noise_fd(default_rng(), model)
function noise_fd(rng::AbstractRNG, model::GaussianNoise)
    out = Matrix{ComplexF64}(undef, length(model.spectrum), 1)
    noise_fd!(rng, out, model)
    return FrequencySeries(out; model.frequency_metadata...)
end
noise_fd!(out::AbstractVecOrMat{<:Complex}, model::GaussianNoise) = noise_fd!(default_rng(), out, model)
function noise_fd!(rng::AbstractRNG, out::AbstractVecOrMat{<:Complex}, model::GaussianNoise)
    size(out, 1) == length(model.spectrum) && size(out, 2) == 1 || throw(DimensionMismatch("noise output shape differs"))
    df = model.frequency_metadata.frequency_spacing
    for i in eachindex(out, model.spectrum)
        scale = sqrt(model.spectrum[i] / df) / 2
        out[i] = scale * complex(randn(rng), randn(rng))
    end
    return out
end
function noise_fd!(rng::AbstractRNG, out::FrequencySeries{<:Complex}, model::GaussianNoise)
    check_compatible(out, model)
    noise_fd!(rng, parent(out), model)
    return out
end
noise_td(model::GaussianNoise) = noise_td(default_rng(), model)
noise_td(rng::AbstractRNG, model::GaussianNoise) = to_timedomain(noise_fd(rng, model); nsamples = model.nsamples)

"""Check stationary sampling/channel compatibility, independently of epoch."""
function check_compatible(s::TimeSeries, model::GaussianNoise)
    nsamples(s) == model.nsamples || throw(DimensionMismatch("noise sample count differs"))
    dt = model.time_metadata.sampling_interval
    check_coordinate(sampling_interval(s), dt, dt, "sample spacings")
    channels(s) == model.time_metadata.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
function check_compatible(s::FrequencySeries, model::GaussianNoise)
    metadata(s).spectrum == :onesided || throw(ArgumentError("noise requires a one-sided spectrum"))
    check_same_frequency_grid(frequencies(s), frequencies(model))
    channels(s) == model.frequency_metadata.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
