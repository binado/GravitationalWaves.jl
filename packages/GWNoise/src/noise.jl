"""
    GaussianNoise(curve::TabulatedPSD, prototype::TimeSeries; kwargs...)
    GaussianNoise(curve::TabulatedPSD; nsamples, epoch, channels, sampling_interval | sampling_frequency,
                  minimum_frequency = first(curve.frequencies), maximum_frequency = nothing)

Prepared finite stationary Gaussian noise segment for one channel.

The one-sided spectrum has `nsamples ÷ 2 + 1` bins. Bins outside
`[minimum_frequency, maximum_frequency]` are zero, as are DC and, for even `nsamples`, the
Nyquist bin. `maximum_frequency = nothing` means the smaller of the end of the PSD
coverage and the Nyquist frequency; the band must lie inside both. Frequency spacing is
`1 / (nsamples * sampling_interval)` and is derived, not stored.

The prototype form takes `nsamples`, `epoch`, the sampling interval and the channel from
the series; keywords override or extend them.
"""
struct GaussianNoise{E <: Real, D <: Real, C <: Tuple}
    spectrum::Vector{Float64}
    nsamples::Int
    epoch::E
    sampling_interval::D
    channels::C
    function GaussianNoise(spectrum::Vector{Float64}, nsamples::Integer, epoch::Real, sampling_interval::Real, channels)
        nsamples >= 2 || throw(ArgumentError("noise needs at least two samples"))
        length(spectrum) == nsamples ÷ 2 + 1 ||
            throw(DimensionMismatch("a spectrum for $nsamples samples has $(nsamples ÷ 2 + 1) bins, got $(length(spectrum))"))
        all(x -> isfinite(x) && x >= 0, spectrum) || throw(ArgumentError("noise spectrum must be finite and nonnegative"))
        spectrum[1] == 0 && (isodd(nsamples) || spectrum[end] == 0) ||
            throw(ArgumentError("the DC and Nyquist bins of a noise spectrum must be zero"))
        tm = GWSeries.time_metadata(; epoch, sampling_interval, channels)
        length(tm.channels) == 1 || throw(ArgumentError("GaussianNoise requires one channel"))
        return new{typeof(tm.epoch), typeof(tm.sampling_interval), typeof(tm.channels)}(
            spectrum, Int(nsamples), tm.epoch, tm.sampling_interval, tm.channels,
        )
    end
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
    nsamples >= 2 || throw(ArgumentError("noise needs at least two samples"))
    fs = inv(tm.sampling_interval)
    fs64 = Float64(fs)
    isfinite(fs64) && fs64 > 0 || throw(ArgumentError("noise sampling frequency must be representable as positive finite Float64"))
    limit = min(last(curve.frequencies), fs / 2)
    upper = something(maximum_frequency, limit)
    first(curve.frequencies) <= minimum_frequency <= upper <= limit ||
        throw(ArgumentError("requested band [$minimum_frequency, $upper] must lie within PSD coverage and Nyquist ($(first(curve.frequencies)), $limit)"))
    df = Float64(inv(nsamples * tm.sampling_interval))
    f = range(0.0; step = df, length = nsamples ÷ 2 + 1)
    spectrum = [minimum_frequency <= x <= upper ? psd(x, curve) : 0.0 for x in f]
    spectrum[1] = 0.0
    iseven(nsamples) && (spectrum[end] = 0.0)
    return GaussianNoise(spectrum, nsamples, tm.epoch, tm.sampling_interval, tm.channels)
end

"""Frequency spacing `1 / (nsamples * sampling_interval)` of the model's spectrum."""
frequency_spacing(model::GaussianNoise) = inv(model.nsamples * model.sampling_interval)
"""Frequency grid of the model's one-sided spectrum, starting at DC."""
frequencies(model::GaussianNoise) = range(0.0; step = Float64(frequency_spacing(model)), length = length(model.spectrum))
"""Copy of the one-sided power spectral density the model samples, zero outside its band."""
psd(model::GaussianNoise) = copy(model.spectrum)

"""
    noise_fd([rng], model) -> FrequencySeries
    noise_fd!([rng], out, model)

Draw a one-channel frequency series with interior variance `PSD / (4Δf)` per component,
using the same normalization as `rfft(strain) * sampling_interval`. `out` may be a vector,
an `nsamples ÷ 2 + 1`-row single-column matrix or a `FrequencySeries` on the model's
frequency grid; a series keeps its own epoch.
"""
noise_fd(model::GaussianNoise) = noise_fd(default_rng(), model)
function noise_fd(rng::AbstractRNG, model::GaussianNoise)
    out = Matrix{ComplexF64}(undef, length(model.spectrum), 1)
    noise_fd!(rng, out, model)
    return FrequencySeries(
        out; epoch = model.epoch, frequency_spacing = frequency_spacing(model),
        channels = model.channels, time_samples = model.nsamples,
    )
end
noise_fd!(out::AbstractVecOrMat{<:Complex}, model::GaussianNoise) = noise_fd!(default_rng(), out, model)
function noise_fd!(rng::AbstractRNG, out::AbstractVecOrMat{<:Complex}, model::GaussianNoise)
    size(out, 1) == length(model.spectrum) && size(out, 2) == 1 ||
        throw(DimensionMismatch("noise output has size $(size(out)), expected ($(length(model.spectrum)), 1)"))
    df = Float64(frequency_spacing(model))
    for i in eachindex(model.spectrum)
        scale = sqrt(model.spectrum[i] / df) / 2
        out[i, 1] = scale * complex(randn(rng), randn(rng))
    end
    return out
end
function noise_fd!(rng::AbstractRNG, out::FrequencySeries{<:Complex}, model::GaussianNoise)
    check_compatible(out, model)
    noise_fd!(rng, parent(out), model)
    return out
end

"""Draw noise as a `TimeSeries`: the inverse transform of [`noise_fd`](@ref). Needs an AbstractFFTs backend such as FFTW."""
noise_td(model::GaussianNoise) = noise_td(default_rng(), model)
noise_td(rng::AbstractRNG, model::GaussianNoise) = to_timedomain(noise_fd(rng, model))

"""Check stationary sampling/channel compatibility, independently of epoch."""
function check_compatible(s::TimeSeries, model::GaussianNoise)
    nsamples(s) == model.nsamples || throw(DimensionMismatch("noise sample count differs"))
    dt = model.sampling_interval
    check_coordinate(sampling_interval(s), dt, dt, "sample spacings")
    channels(s) == model.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
function check_compatible(s::FrequencySeries, model::GaussianNoise)
    check_same_frequency_grid(frequencies(s), frequencies(model))
    known = time_samples(s)
    known === nothing || known == model.nsamples || throw(ArgumentError("noise time length differs: $known vs $(model.nsamples)"))
    channels(s) == model.channels || throw(ArgumentError("noise channel orders differ"))
    return nothing
end
