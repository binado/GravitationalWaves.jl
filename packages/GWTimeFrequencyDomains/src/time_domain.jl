"""TimeDomain(; epoch, nsamples, sampling_frequency, channels): immutable sampling metadata."""
struct TimeDomain{E <: Real, D <: Real, C <: Tuple} <: AbstractDomain
    epoch::E
    nsamples::Int
    spacing::D
    channels::C
    TimeDomain(epoch, nsamples, spacing, channels, ::Validated) =
        new{typeof(epoch), typeof(spacing), typeof(channels)}(epoch, nsamples, spacing, channels)
end

function TimeDomain(; epoch::Real, nsamples::Integer, sampling_frequency::Real, channels)
    t0, fs = metadata(epoch), metadata(sampling_frequency)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    nsamples >= 0 || throw(ArgumentError("sample count must be nonnegative"))
    isfinite(fs) && fs > 0 || throw(ArgumentError("sampling frequency must be finite and positive"))
    dt = inv(fs)
    isfinite(dt) && dt > 0 || throw(ArgumentError("sample spacing must be finite and positive"))
    return TimeDomain(t0, Int(nsamples), dt, channel_order(channels), Validated())
end

"""Stored time sample count."""
nsamples(d::TimeDomain) = d.nsamples
"""Time spacing in seconds."""
sample_spacing(d::TimeDomain) = d.spacing
"""Sampling frequency in Hz, derived from time spacing."""
sampling_frequency(d::TimeDomain) = inv(sample_spacing(d))
"""Half-open time segment duration N*dt."""
duration(d::TimeDomain) = nsamples(d) * sample_spacing(d)
"""Lazy absolute sample times."""
times(d::TimeDomain) = range(epoch(d); step = sample_spacing(d), length = nsamples(d))

function check_fourier_domain(d::TimeDomain)
    nsamples(d) >= 2 && iseven(nsamples(d)) || throw(ArgumentError("forward FFT needs even N >= 2"))
    return nothing
end

"""Derive the one-sided Fourier domain; samples are not transformed."""
function frequency_domain(d::TimeDomain)
    check_fourier_domain(d)
    f = range(0.0; step = sampling_frequency(d) / nsamples(d), length = nsamples(d) ÷ 2 + 1)
    return FrequencyDomain(; frequencies = f, epoch = epoch(d), channels = channels(d))
end

coordinates(d::TimeDomain) = times(d)
rebuild(d::TimeDomain, n, origin) = TimeDomain(origin, n, sample_spacing(d), channels(d), Validated())
