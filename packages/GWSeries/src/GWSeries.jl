"""Explicit sampling containers with array kernel delegation."""
module GWSeries

export TimeSeries, FrequencySeries, epoch, channels, nsamples, times, frequencies
export sample_spacing, sampling_frequency, duration, frequency_spacing
export with_data, crop, pad, check_same_grid, to_frequencyseries, to_timeseries

abstract type AbstractSeries end
const Samples = Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}

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
struct TimeSeries{A, E, D, C} <: AbstractSeries
    data::A
    epoch::E
    spacing::D
    channels::C
    TimeSeries(data, epoch, spacing, channels, ::Val{:validated}) =
        new{typeof(data), typeof(epoch), typeof(spacing), typeof(channels)}(data, epoch, spacing, channels)
end
function TimeSeries(data::Samples; epoch::Real, sampling_frequency::Real, channels)
    t0, fs = metadata(epoch), metadata(sampling_frequency)
    names = channel_order(channels)
    validate(data, t0, names)
    isfinite(fs) && fs > 0 || throw(ArgumentError("sampling frequency must be finite and positive"))
    dt = inv(fs)
    isfinite(dt) && dt > 0 || throw(ArgumentError("sample spacing must be finite and positive"))
    return TimeSeries(data, t0, dt, names, Val(:validated))
end

"""FrequencySeries(data; frequencies, epoch, channels): uniformly spaced, nonnegative bins."""
struct FrequencySeries{A, E, F, C} <: AbstractSeries
    data::A
    epoch::E
    frequencies::F
    channels::C
    FrequencySeries(data, epoch, frequencies, channels, ::Val{:validated}) =
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
    return FrequencySeries(data, t0, f, names, Val(:validated))
end

Base.parent(x::AbstractSeries) = x.data
Base.copy(x::AbstractSeries) = with_data(x, copy(parent(x)))
epoch(x::AbstractSeries) = x.epoch
channels(x::AbstractSeries) = x.channels
nsamples(x::AbstractSeries) = size(parent(x), 1)
sample_spacing(x::TimeSeries) = x.spacing
sampling_frequency(x::TimeSeries) = inv(sample_spacing(x))
duration(x::TimeSeries) = nsamples(x) * sample_spacing(x)
times(x::TimeSeries) = range(epoch(x); step = sample_spacing(x), length = nsamples(x))
frequencies(x::FrequencySeries) = x.frequencies
frequency_spacing(x::FrequencySeries) = step(frequencies(x))
coordinates(x::TimeSeries) = times(x)
coordinates(x::FrequencySeries) = frequencies(x)
rebuild(x::TimeSeries, data, origin) = TimeSeries(data, origin, sample_spacing(x), channels(x), Val(:validated))
rebuild(x::FrequencySeries, data, origin) = FrequencySeries(data, epoch(x), range(origin; step = frequency_spacing(x), length = size(data, 1)), channels(x), Val(:validated))

"""Replace storage without copying; require identical complete shape and preserve metadata."""
function with_data(x::AbstractSeries, data::Samples)
    Base.require_one_based_indexing(data)
    size(data) == size(parent(x)) || throw(DimensionMismatch("replacement data must have identical shape"))
    return rebuild(x, data, first(coordinates(x)))
end

# Four ulps of each absolute coordinate and eight ulps of relative displacement/spacing.
ulp(x) = eps(float(x))
function coordinate_tolerance(a, b, spacing)
    tol = 4 * (ulp(a) + ulp(b)) + 8 * (ulp(a - b) + ulp(spacing))
    tol < spacing / 2 || throw(ArgumentError("coordinate precision cannot resolve this grid"))
    return tol
end
function check_coordinate(a, b, spacing, label)
    abs(a - b) <= coordinate_tolerance(a, b, spacing) || throw(ArgumentError("$label differ"))
    return nothing
end
"""Validate equal domain, grid, epoch, and ordered channels; return nothing."""
function check_same_grid(a::AbstractSeries, b::AbstractSeries)
    (a isa TimeSeries) == (b isa TimeSeries) || throw(ArgumentError("series domains differ"))
    nsamples(a) == nsamples(b) || throw(DimensionMismatch("sample counts differ"))
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    fa, fb = coordinates(a), coordinates(b)
    check_frequency_grid(fa, fb)
    epoch_spacing = a isa TimeSeries ? step(fa) : inv(max(2 * (nsamples(a) - 1), 1) * step(fa))
    check_coordinate(epoch(a), epoch(b), epoch_spacing, "epochs")
    return nothing
end
# Shared by stationary consumers, which deliberately have no epoch.
function check_frequency_grid(a::AbstractRange, b::AbstractRange)
    length(a) == length(b) || throw(DimensionMismatch("grid sample counts differ"))
    check_coordinate(step(a), step(b), min(step(a), step(b)), "spacings")
    check_coordinate(first(a), first(b), step(a), "origins")
    check_coordinate(length(a) * step(a), length(b) * step(b), step(a), "grid extents")
    return nothing
end

function boundary_index(x, boundary::Real)
    grid = coordinates(x)
    origin, spacing = first(grid), step(grid)
    value = metadata(boundary)
    isfinite(value) || throw(ArgumentError("crop boundary must be finite"))
    displacement = value - origin
    tol = coordinate_tolerance(value, origin, spacing)
    -tol <= displacement <= nsamples(x) * spacing + tol || throw(ArgumentError("crop boundary outside grid extent"))
    index = round(Int, displacement / spacing)
    abs(displacement - index * spacing) <= tol || throw(ArgumentError("crop boundary is off grid"))
    return clamp(index, 0, nsamples(x))
end
"""Copy a half-open interval with grid-aligned boundaries."""
function crop(x::AbstractSeries; start::Real, stop::Real)
    lo, hi = boundary_index(x, start), boundary_index(x, stop)
    lo <= hi || throw(ArgumentError("crop start must precede stop"))
    data = ndims(parent(x)) == 1 ? parent(x)[(lo + 1):hi] : parent(x)[(lo + 1):hi, :]
    return rebuild(x, data, first(coordinates(x)) + lo * step(coordinates(x)))
end
"""Zero pad by sample/bin counts; shift the origin when padding before."""
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

"""Physical dt*rfft conversion, even N ≥ 2; load FFTW for CPU implementations."""
function to_frequencyseries end
"""Inverse physical FFT on a complete zero-origin grid; imaginary endpoints are discarded."""
function to_timeseries end

end
