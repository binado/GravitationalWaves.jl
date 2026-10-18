"""
    GaussianNoise(curve; sampling_frequency, nsamples,
                  minimum_frequency=first(curve.frequencies),
                  maximum_frequency=min(last(curve.frequencies), sampling_frequency/2))

Prepare a finite stationary Gaussian noise segment. `nsamples` and
`sampling_frequency` (Hz) determine the duration and one-sided FFT grid.
The active band must be covered by the curve; other bins have zero power.
DC and, for even sample counts, Nyquist are zero. No taper is applied.
Draws are independent periodic FFT segments, not a continuous noise stream.
"""
struct GaussianNoise{F <: AbstractRange}
    sampling_frequency::Float64
    nsamples::Int
    frequencies::F
    spectrum::Vector{Float64}
    scales::Vector{Float64}
end

function GaussianNoise(
        curve::TabulatedPSD; sampling_frequency, nsamples::Integer,
        minimum_frequency = first(curve.frequencies),
        maximum_frequency = min(last(curve.frequencies), sampling_frequency / 2),
    )
    isfinite(sampling_frequency) && sampling_frequency > 0 && nsamples >= 2 ||
        throw(ArgumentError("need a positive sampling frequency and at least two samples"))
    first(curve.frequencies) <= minimum_frequency <= maximum_frequency <=
        min(last(curve.frequencies), sampling_frequency / 2) ||
        throw(ArgumentError("requested band must lie within PSD coverage and Nyquist"))
    fs = Float64(sampling_frequency)
    df = fs / nsamples
    f = range(0.0; step = df, length = nsamples ÷ 2 + 1)
    s = [minimum_frequency <= x <= maximum_frequency ? psd(x, curve) : 0.0 for x in f]
    s[1] = 0.0
    iseven(nsamples) && (s[end] = 0.0)
    scales = sqrt.(s ./ df) ./ 2
    return GaussianNoise(fs, Int(nsamples), f, s, scales)
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
    destination = Vector{ComplexF64}(undef, length(model.frequencies))
    noise_fd!(rng, destination, model)
    return destination
end

"""`noise_fd!([rng], destination, model)` fills a complex vector with a new FD draw."""
noise_fd!(destination::AbstractVector{<:Complex}, model::GaussianNoise) =
    noise_fd!(default_rng(), destination, model)

function noise_fd!(rng::AbstractRNG, destination::AbstractVector{<:Complex}, model::GaussianNoise)
    length(destination) == length(model.scales) || throw(DimensionMismatch("wrong FD output length"))
    for i in eachindex(destination, model.scales)
        destination[i] = model.scales[i] * complex(randn(rng), randn(rng))
    end
    return destination
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
    noise_td!([rng], destination, model)

Fill a real vector with a new TD draw. This convenience method allocates FFT
temporaries; it does not cache scratch buffers or FFT plans in the model.
"""
noise_td!(destination::AbstractVector{<:Real}, model::GaussianNoise) =
    noise_td!(default_rng(), destination, model)

function noise_td!(rng::AbstractRNG, destination::AbstractVector{<:Real}, model::GaussianNoise)
    length(destination) == model.nsamples || throw(DimensionMismatch("wrong TD output length"))
    return copyto!(destination, noise_td(rng, model))
end

"""Draw single-channel series using an explicit epoch and channel label."""
noise_fd(::Type{FrequencySeries}, model::GaussianNoise; kwargs...) =
    noise_fd(default_rng(), FrequencySeries, model; kwargs...)
function noise_fd(rng::AbstractRNG, ::Type{FrequencySeries}, model::GaussianNoise; epoch, channels)
    iseven(model.nsamples) || throw(ArgumentError("FD series noise requires an even-length noise model"))
    # Validate metadata before drawing, including channel count.
    prototype = FrequencySeries(model.scales; frequencies = model.frequencies, epoch, channels)
    return with_data(prototype, noise_fd(rng, model))
end

"""Draw single-channel time series; odd-length segments remain supported."""
noise_td(::Type{TimeSeries}, model::GaussianNoise; kwargs...) =
    noise_td(default_rng(), TimeSeries, model; kwargs...)
function noise_td(rng::AbstractRNG, ::Type{TimeSeries}, model::GaussianNoise; epoch, channels)
    # A lazy range validates metadata/shape without allocating another sample buffer.
    prototype = TimeSeries(1:model.nsamples; sampling_frequency = model.sampling_frequency, epoch, channels)
    return with_data(prototype, noise_td(rng, model))
end

noise_fd!(destination::FrequencySeries, model::GaussianNoise) = noise_fd!(default_rng(), destination, model)
function noise_fd!(rng::AbstractRNG, destination::FrequencySeries, model::GaussianNoise)
    parent(destination) isa AbstractVector{<:Complex} || throw(ArgumentError("FD noise destination must hold a complex vector"))
    iseven(model.nsamples) || throw(ArgumentError("FD series noise requires an even-length noise model"))
    prototype = FrequencySeries(parent(destination); frequencies = model.frequencies, epoch = epoch(destination), channels = channels(destination))
    check_same_grid(destination, prototype)
    noise_fd!(rng, parent(destination), model)
    return destination
end

noise_td!(destination::TimeSeries, model::GaussianNoise) = noise_td!(default_rng(), destination, model)
function noise_td!(rng::AbstractRNG, destination::TimeSeries, model::GaussianNoise)
    parent(destination) isa AbstractVector{<:Real} || throw(ArgumentError("TD noise destination must hold a real vector"))
    length(parent(destination)) == model.nsamples || throw(DimensionMismatch("wrong TD output length"))
    prototype = TimeSeries(parent(destination); sampling_frequency = model.sampling_frequency, epoch = epoch(destination), channels = channels(destination))
    check_same_grid(destination, prototype)
    noise_td!(rng, parent(destination), model)
    return destination
end
