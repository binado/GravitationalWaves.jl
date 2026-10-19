"""FrequencyDomain(; frequencies, epoch, channels): uniform nonnegative frequency grid metadata."""
struct FrequencyDomain{E <: Real, F <: AbstractRange, C <: Tuple} <: AbstractDomain
    epoch::E
    frequencies::F
    channels::C
    FrequencyDomain(epoch, frequencies, channels, ::Validated) =
        new{typeof(epoch), typeof(frequencies), typeof(channels)}(epoch, frequencies, channels)
end

function FrequencyDomain(; frequencies::AbstractRange{<:Real}, epoch::Real, channels)
    t0 = metadata(epoch)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    f0, df = metadata(first(frequencies)), metadata(step(frequencies))
    isfinite(f0) && f0 >= 0 && isfinite(df) && df > 0 ||
        throw(ArgumentError("frequencies need a finite nonnegative origin and positive spacing"))
    f = range(f0; step = df, length = length(frequencies))
    isempty(f) || isfinite(last(f)) || throw(ArgumentError("frequencies must be finite"))
    return FrequencyDomain(t0, f, channel_order(channels), Validated())
end

"""Stored frequency bin count."""
nsamples(d::FrequencyDomain) = length(d.frequencies)
"""Lazy bin frequencies."""
frequencies(d::FrequencyDomain) = d.frequencies
"""Frequency spacing in Hz."""
frequency_spacing(d::FrequencyDomain) = step(frequencies(d))

function check_fourier_domain(d::FrequencyDomain)
    nsamples(d) >= 2 && first(frequencies(d)) == 0 ||
        throw(ArgumentError("inverse FFT needs at least two bins starting at zero; restore omitted low bins explicitly"))
    return nothing
end

"""Infer the even-length inverse time domain N=2(M-1)."""
function time_domain(d::FrequencyDomain)
    check_fourier_domain(d)
    n = 2 * (nsamples(d) - 1)
    return TimeDomain(; epoch = epoch(d), nsamples = n, sampling_frequency = n * frequency_spacing(d), channels = channels(d))
end

coordinates(d::FrequencyDomain) = frequencies(d)
rebuild(d::FrequencyDomain, n, origin) = FrequencyDomain(epoch(d), range(origin; step = frequency_spacing(d), length = n), channels(d), Validated())
