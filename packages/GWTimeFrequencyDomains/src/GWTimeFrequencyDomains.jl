"""Detached sampling domains, explicit validation, and physical Fourier conversion."""
module GWTimeFrequencyDomains

export TimeDomain, FrequencyDomain, epoch, channels, nsamples, times, frequencies
export sample_spacing, sampling_frequency, duration, frequency_spacing
export check_data, check_same_grid, check_same_frequency_grid, check_compatible, check_fourier_domain
export crop, pad, frequency_domain, time_domain, to_frequencydomain, to_timedomain
export plan_to_frequencydomain, plan_to_timedomain, output_domain, transform, transform!

abstract type AbstractDomain end
const Samples = Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}
struct Validated end
metadata(x::Real) = float(promote(x, 0.0)[1])
function channel_order(names)
    labels = Tuple(names)
    !isempty(labels) && all(x -> x isa Symbol || x isa AbstractString, labels) ||
        throw(ArgumentError("channels must be nonempty symbol or string names"))
    canonical = map(Symbol, labels)
    allunique(canonical) || throw(ArgumentError("channel names must be unique"))
    return canonical
end

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

"""Segment epoch in seconds, independent of sample precision."""
epoch(d::AbstractDomain) = d.epoch
"""Ordered, unique Symbol channel names."""
channels(d::AbstractDomain) = d.channels
"""Stored time sample count or frequency bin count."""
nsamples(d::TimeDomain) = d.nsamples
nsamples(d::FrequencyDomain) = length(d.frequencies)
"""Time spacing in seconds."""
sample_spacing(d::TimeDomain) = d.spacing
"""Sampling frequency in Hz, derived from time spacing."""
sampling_frequency(d::TimeDomain) = inv(sample_spacing(d))
"""Half-open time segment duration N*dt."""
duration(d::TimeDomain) = nsamples(d) * sample_spacing(d)
"""Lazy absolute sample times."""
times(d::TimeDomain) = range(epoch(d); step = sample_spacing(d), length = nsamples(d))
"""Lazy bin frequencies."""
frequencies(d::FrequencyDomain) = d.frequencies
"""Frequency spacing in Hz."""
frequency_spacing(d::FrequencyDomain) = step(frequencies(d))
coordinates(d::TimeDomain) = times(d)
coordinates(d::FrequencyDomain) = frequencies(d)
rebuild(d::TimeDomain, n, origin) = TimeDomain(origin, n, sample_spacing(d), channels(d), Validated())
rebuild(d::FrequencyDomain, n, origin) = FrequencyDomain(epoch(d), range(origin; step = frequency_spacing(d), length = n), channels(d), Validated())

"""Check numeric vector/matrix storage, one-based axes, sample count and channel count on the host."""
function check_data(d::AbstractDomain, values::Samples)
    Base.require_one_based_indexing(values)
    size(values, 1) == nsamples(d) || throw(DimensionMismatch("sample count does not match domain"))
    size(values, 2) == length(channels(d)) || throw(DimensionMismatch("channel count does not match domain"))
    return nothing
end
check_data(::AbstractDomain, values) = throw(ArgumentError("data must be a numeric vector or matrix, got $(typeof(values))"))

"""Host compatibility check extended by packages owning prepared numerical objects."""
function check_compatible end

"""Check even-length one-sided real-FFT eligibility on the host."""
function check_fourier_domain(d::TimeDomain)
    nsamples(d) >= 2 && iseven(nsamples(d)) || throw(ArgumentError("forward FFT needs even N >= 2"))
    return nothing
end
function check_fourier_domain(d::FrequencyDomain)
    nsamples(d) >= 2 && first(frequencies(d)) == 0 ||
        throw(ArgumentError("inverse FFT needs at least two bins starting at zero; restore omitted low bins explicitly"))
    return nothing
end
"""Derive and validate the one-sided Fourier domain; samples are not transformed."""
function frequency_domain(d::TimeDomain)
    check_fourier_domain(d)
    f = range(0.0; step = sampling_frequency(d) / nsamples(d), length = nsamples(d) ÷ 2 + 1)
    return FrequencyDomain(; frequencies = f, epoch = epoch(d), channels = channels(d))
end
"""Infer and validate the even-length inverse time domain N=2(M-1)."""
function time_domain(d::FrequencyDomain)
    check_fourier_domain(d)
    n = 2 * (nsamples(d) - 1)
    return TimeDomain(; epoch = epoch(d), nsamples = n, sampling_frequency = n * frequency_spacing(d), channels = channels(d))
end
# Four ulps of each absolute coordinate and eight ulps of relative displacement/spacing.
# At GPS magnitude (1.2e9 s) this is about 2 microseconds for epochs; Float64 cannot
# resolve finer offsets, so use BigFloat metadata when sub-microsecond timing matters.
ulp(x) = eps(float(x))
function coordinate_tolerance(a, b, spacing, label)
    tol = 4 * (ulp(a) + ulp(b)) + 8 * (ulp(a - b) + ulp(spacing))
    tol < spacing / 2 ||
        throw(ArgumentError("$label: coordinate precision cannot resolve this grid (tolerance $tol, spacing $spacing)"))
    return tol
end
function check_coordinate(a, b, spacing, label, hint = "")
    a == b && return nothing
    tol = coordinate_tolerance(a, b, spacing, label)
    abs(a - b) <= tol || throw(ArgumentError("$label differ: $a vs $b (tolerance $tol)$hint"))
    return nothing
end
const EPOCH_HINT = "; frequency-domain phases are referenced to the epoch, so build both domains with the same start time"

# Count, spacing, origin, and accumulated extent; `origin` names the origin in errors.
function check_range(a::AbstractRange, b::AbstractRange, origin, hint = "")
    length(a) == length(b) || throw(DimensionMismatch("grid sample counts differ"))
    check_coordinate(step(a), step(b), min(step(a), step(b)), "spacings")
    check_coordinate(first(a), first(b), step(a), origin, hint)
    check_coordinate(length(a) * step(a), length(b) * step(b), step(a), "grid extents")
    return nothing
end

"""
    check_same_grid(a, b)

Validate equal domain, sample count, grid, epoch, and ordered channels; return nothing.
Throws `ArgumentError` (or `DimensionMismatch` for counts) naming the differing
quantity, both values, and the tolerance.

Coordinates are compared with a few-ulp tolerance (see the README), not exactly;
identical values always pass. Time domains compare the epoch as the grid origin.
Frequency domains additionally compare epochs using the spacing `1/(N df)` of the
even-length FFT that would produce them. Epochs must match because
frequency-domain phases are referenced to the epoch.
"""
function check_same_grid(a::TimeDomain, b::TimeDomain)
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    check_range(times(a), times(b), "epochs", EPOCH_HINT)
    return nothing
end
function check_same_grid(a::FrequencyDomain, b::FrequencyDomain)
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    check_range(frequencies(a), frequencies(b), "origins")
    epoch_spacing = inv(max(2 * (nsamples(a) - 1), 1) * frequency_spacing(a))
    check_coordinate(epoch(a), epoch(b), epoch_spacing, "epochs", EPOCH_HINT)
    return nothing
end
check_same_grid(::AbstractDomain, ::AbstractDomain) = throw(ArgumentError("domain kinds differ"))

"""
    check_same_frequency_grid(a, b)

Check that two frequency ranges have equal counts, spacings, origins, and extents
within the [`check_same_grid`](@ref) tolerance; return nothing. For stationary
consumers such as noise metrics, which deliberately have no epoch.
"""
check_same_frequency_grid(a::AbstractRange, b::AbstractRange) = check_range(a, b, "origins")

function boundary_index(x, boundary::Real)
    grid = coordinates(x)
    origin, spacing = first(grid), step(grid)
    value = metadata(boundary)
    isfinite(value) || throw(ArgumentError("crop boundary must be finite"))
    displacement = value - origin
    value == origin && return 0
    tol = coordinate_tolerance(value, origin, spacing, "crop boundary")
    -tol <= displacement <= nsamples(x) * spacing + tol || throw(ArgumentError("crop boundary outside grid extent"))
    index = round(Int, displacement / spacing)
    abs(displacement - index * spacing) <= tol || throw(ArgumentError("crop boundary is off grid"))
    return clamp(index, 0, nsamples(x))
end
"""
    crop(domain, values; start, stop)

Copy the half-open interval `[start, stop)` with grid-aligned boundaries. Bounds
are in seconds for a `TimeDomain` and Hz for a `FrequencyDomain`. Throws
`ArgumentError` for non-finite, off-grid, or out-of-extent boundaries and when
`start` follows `stop`.
"""
function crop(x::AbstractDomain, values::Samples; start::Real, stop::Real)
    check_data(x, values)
    lo, hi = boundary_index(x, start), boundary_index(x, stop)
    lo <= hi || throw(ArgumentError("crop start must precede stop"))
    data = ndims(values) == 1 ? values[(lo + 1):hi] : values[(lo + 1):hi, :]
    return (rebuild(x, size(data, 1), first(coordinates(x)) + lo * step(coordinates(x))), data)
end
"""
    pad(domain, values; before = 0, after = 0)

Zero pad by sample/bin counts; shift the epoch (time) or origin (frequency) when
padding before. Counts must be nonnegative, and frequency padding cannot create
negative frequencies.
"""
function pad(x::AbstractDomain, values::Samples; before::Integer = 0, after::Integer = 0)
    check_data(x, values)
    before >= 0 && after >= 0 || throw(ArgumentError("padding counts must be nonnegative"))
    n = Base.checked_add(Base.checked_add(nsamples(x), Int(before)), Int(after))
    origin = first(coordinates(x)) - before * step(coordinates(x))
    isfinite(origin) || throw(ArgumentError("padded origin must be finite"))
    x isa FrequencyDomain && origin < 0 && throw(ArgumentError("padding creates negative frequencies"))
    isfinite(origin + max(n - 1, 0) * step(coordinates(x))) ||
        throw(ArgumentError("padded coordinates must be finite"))
    dims = ndims(values) == 1 ? (n,) : (n, size(values, 2))
    data = similar(values, dims)
    fill!(data, zero(eltype(data)))
    if ndims(data) == 1
        data[(before + 1):(before + nsamples(x))] .= values
    else
        data[(before + 1):(before + nsamples(x)), :] .= values
    end
    return (rebuild(x, size(data, 1), origin), data)
end


"""Physical dt*rfft of validated real CPU samples along dimension one. Load FFTW. No domain checks."""
function to_frequencydomain end
"""Physical inverse FFT of validated CPU bins. Imaginary endpoints are discarded. No domain checks."""
function to_timedomain end
"""Prepare a forward FFT on the host, validating the domain and prototype."""
function plan_to_frequencydomain end
"""Prepare an inverse FFT on the host, validating the domain and prototype."""
function plan_to_timedomain end
"""The prepared transform's immutable output domain."""
function output_domain end
"""Allocate and execute a prepared physical transform without domain validation."""
function transform end
"""Execute into separate destination storage; inputs must match the prepared prototype."""
function transform! end

end
