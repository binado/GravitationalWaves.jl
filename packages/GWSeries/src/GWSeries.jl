"""
Explicit sampling containers with array kernel delegation.

Series wrap an array with its grid, epoch, and channel labels. They do not
subtype `AbstractArray`; use the accessors below (`parent`, `epoch`, `channels`,
...) rather than struct fields, and unwrap to arrays before numerical kernels.
`parent`, `with_data`, and the constructors alias their storage, `copy` and
`crop` copy. Mutating the aliased array is visible through the series; resizing
it breaks the shape invariant and is unsupported.
"""
module GWSeries

export TimeSeries, FrequencySeries, epoch, channels, nsamples, times, frequencies
export sample_spacing, sampling_frequency, duration, frequency_spacing
export with_data, crop, pad, check_same_grid, check_same_frequency_grid
export to_frequencyseries, to_timeseries

abstract type AbstractSeries end
const Samples = Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}

# Marks constructor calls whose arguments were already validated. Internal: call
# the keyword constructors instead.
struct Validated end

# Metadata precision is independent of sample precision.
metadata(x::Real) = float(promote(x, 0.0)[1])
function channel_order(names)
    labels = Tuple(names)
    !isempty(labels) && all(x -> x isa Symbol || x isa AbstractString, labels) ||
        throw(ArgumentError("channels must be nonempty symbol or string names"))
    canonical = map(Symbol, labels)
    allunique(canonical) || throw(ArgumentError("channel names must be unique"))
    return canonical
end
function validate(data, t0, names)
    Base.require_one_based_indexing(data)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    size(data, 2) == length(names) || throw(DimensionMismatch("channels must match sample columns"))
    return nothing
end

"""TimeSeries(data; epoch, sampling_frequency, channels): retain numeric storage without copying."""
struct TimeSeries{A <: Samples, E <: Real, D <: Real, C <: Tuple} <: AbstractSeries
    data::A
    epoch::E
    spacing::D
    channels::C
    TimeSeries(data, epoch, spacing, channels, ::Validated) =
        new{typeof(data), typeof(epoch), typeof(spacing), typeof(channels)}(data, epoch, spacing, channels)
end
function TimeSeries(data::Samples; epoch::Real, sampling_frequency::Real, channels)
    t0, fs = metadata(epoch), metadata(sampling_frequency)
    names = channel_order(channels)
    validate(data, t0, names)
    isfinite(fs) && fs > 0 || throw(ArgumentError("sampling frequency must be finite and positive"))
    dt = inv(fs)
    isfinite(dt) && dt > 0 || throw(ArgumentError("sample spacing must be finite and positive"))
    return TimeSeries(data, t0, dt, names, Validated())
end

"""FrequencySeries(data; frequencies, epoch, channels): uniformly spaced, nonnegative bins."""
struct FrequencySeries{A <: Samples, E <: Real, F <: AbstractRange, C <: Tuple} <: AbstractSeries
    data::A
    epoch::E
    frequencies::F
    channels::C
    FrequencySeries(data, epoch, frequencies, channels, ::Validated) =
        new{typeof(data), typeof(epoch), typeof(frequencies), typeof(channels)}(data, epoch, frequencies, channels)
end
function FrequencySeries(data::Samples; frequencies::AbstractRange{<:Real}, epoch::Real, channels)
    t0 = metadata(epoch)
    names = channel_order(channels)
    validate(data, t0, names)
    f0, df = metadata(first(frequencies)), metadata(step(frequencies))
    isfinite(f0) && f0 >= 0 && isfinite(df) && df > 0 ||
        throw(ArgumentError("frequencies need a finite nonnegative origin and positive spacing"))
    size(data, 1) == length(frequencies) || throw(DimensionMismatch("frequency count must match samples"))
    f = range(f0; step = df, length = length(frequencies))
    isempty(f) || isfinite(last(f)) || throw(ArgumentError("frequencies must be finite"))
    return FrequencySeries(data, t0, f, names, Validated())
end

"""The wrapped array itself (not a copy); resizing it is unsupported."""
Base.parent(x::AbstractSeries) = x.data
"""Copy the samples; the metadata is immutable and shared."""
Base.copy(x::AbstractSeries) = with_data(x, copy(parent(x)))
"""Epoch in seconds (GPS for detector data); the segment start for frequency series."""
epoch(x::AbstractSeries) = x.epoch
"""Canonical `Symbol` channel labels, one per sample column."""
channels(x::AbstractSeries) = x.channels
"""Number of rows (time samples or frequency bins)."""
nsamples(x::AbstractSeries) = size(parent(x), 1)
"""Time between samples, `inv(sampling_frequency)`, stored as given at construction."""
sample_spacing(x::TimeSeries) = x.spacing
"""Sampling frequency, `inv(sample_spacing(x))`; may differ from the input by an ulp."""
sampling_frequency(x::TimeSeries) = inv(sample_spacing(x))
"""Half-open segment length `nsamples * sample_spacing` of a time series."""
duration(x::TimeSeries) = nsamples(x) * sample_spacing(x)
"""Lazy range of sample times starting at the epoch."""
times(x::TimeSeries) = range(epoch(x); step = sample_spacing(x), length = nsamples(x))
"""Lazy `Float64`-or-wider range of bin frequencies in Hz."""
frequencies(x::FrequencySeries) = x.frequencies
"""Bin width in Hz."""
frequency_spacing(x::FrequencySeries) = step(frequencies(x))
coordinates(x::TimeSeries) = times(x)
coordinates(x::FrequencySeries) = frequencies(x)
rebuild(x::TimeSeries, data, origin) = TimeSeries(data, origin, sample_spacing(x), channels(x), Validated())
rebuild(x::FrequencySeries, data, origin) = FrequencySeries(data, epoch(x), range(origin; step = frequency_spacing(x), length = size(data, 1)), channels(x), Validated())

"""Replace storage without copying; require identical complete shape and preserve metadata."""
function with_data(x::AbstractSeries, data::Samples)
    Base.require_one_based_indexing(data)
    size(data) == size(parent(x)) || throw(DimensionMismatch("replacement data must have identical shape"))
    return rebuild(x, data, first(coordinates(x)))
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
const EPOCH_HINT = "; frequency-domain phases are referenced to the epoch, so build both series with the same start time"

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
identical values always pass. Time series compare the epoch as the grid origin.
Frequency series additionally compare epochs using the spacing `1/(N df)` of the
even-length FFT that would produce them. Epochs must match because
frequency-domain phases are referenced to the epoch.
"""
function check_same_grid(a::AbstractSeries, b::AbstractSeries)
    (a isa TimeSeries) == (b isa TimeSeries) || throw(ArgumentError("series domains differ"))
    nsamples(a) == nsamples(b) || throw(DimensionMismatch("sample counts differ"))
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    if a isa TimeSeries
        check_range(times(a), times(b), "epochs", EPOCH_HINT)
    else
        check_range(frequencies(a), frequencies(b), "origins")
        epoch_spacing = inv(max(2 * (nsamples(a) - 1), 1) * frequency_spacing(a))
        check_coordinate(epoch(a), epoch(b), epoch_spacing, "epochs", EPOCH_HINT)
    end
    return nothing
end
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
    crop(x; start, stop)

Copy the half-open interval `[start, stop)` with grid-aligned boundaries. Bounds
are in seconds for a `TimeSeries` and Hz for a `FrequencySeries`. Throws
`ArgumentError` for non-finite, off-grid, or out-of-extent boundaries and when
`start` follows `stop`.
"""
function crop(x::AbstractSeries; start::Real, stop::Real)
    lo, hi = boundary_index(x, start), boundary_index(x, stop)
    lo <= hi || throw(ArgumentError("crop start must precede stop"))
    data = ndims(parent(x)) == 1 ? parent(x)[(lo + 1):hi] : parent(x)[(lo + 1):hi, :]
    return rebuild(x, data, first(coordinates(x)) + lo * step(coordinates(x)))
end
"""
    pad(x; before = 0, after = 0)

Zero pad by sample/bin counts; shift the epoch (time) or origin (frequency) when
padding before. Counts must be nonnegative, and frequency padding cannot create
negative frequencies.
"""
function pad(x::AbstractSeries; before::Integer = 0, after::Integer = 0)
    before >= 0 && after >= 0 || throw(ArgumentError("padding counts must be nonnegative"))
    n = Base.checked_add(Base.checked_add(nsamples(x), Int(before)), Int(after))
    origin = first(coordinates(x)) - before * step(coordinates(x))
    isfinite(origin) || throw(ArgumentError("padded origin must be finite"))
    x isa FrequencySeries && origin < 0 && throw(ArgumentError("padding creates negative frequencies"))
    isfinite(origin + max(n - 1, 0) * step(coordinates(x))) ||
        throw(ArgumentError("padded coordinates must be finite"))
    dims = ndims(parent(x)) == 1 ? (n,) : (n, size(parent(x), 2))
    data = similar(parent(x), dims)
    fill!(data, zero(eltype(data)))
    if ndims(data) == 1
        data[(before + 1):(before + nsamples(x))] .= parent(x)
    else
        data[(before + 1):(before + nsamples(x)), :] .= parent(x)
    end
    return rebuild(x, data, origin)
end

"""
    to_frequencyseries(x::TimeSeries)

Physical `dt * rfft` conversion of real CPU data, even N ≥ 2. Load `FFTW` to
activate it; other inputs throw `ArgumentError`.
"""
to_frequencyseries(x::AbstractSeries) =
    unsupported_conversion(:to_frequencyseries, "a TimeSeries of Float32/Float64 CPU array storage", x)
"""
    to_timeseries(h::FrequencySeries)

Inverse physical FFT of a CPU spectrum on a complete zero-origin grid.
The imaginary parts of the DC and Nyquist bins are discarded. Load `FFTW` to
activate it; other inputs throw `ArgumentError`.
"""
to_timeseries(h::AbstractSeries) =
    unsupported_conversion(:to_timeseries, "a FrequencySeries of Float32/Float64/ComplexF32/ComplexF64 CPU array storage", h)

function unsupported_conversion(name, expected, x)
    Base.get_extension(GWSeries, :GWSeriesFFTWExt) === nothing &&
        throw(ArgumentError("$name needs the FFTW extension; run `using FFTW`"))
    throw(ArgumentError("$name needs $expected (arrays or strided views), got $(typeof(x))"))
end

end
