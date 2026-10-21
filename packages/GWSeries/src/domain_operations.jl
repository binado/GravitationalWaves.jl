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
function check_same_grid(a::TimeSeries, b::TimeSeries)
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    check_range(times(a), times(b), "epochs", EPOCH_HINT)
    return nothing
end
function check_same_grid(a::FrequencySeries, b::FrequencySeries)
    channels(a) == channels(b) || throw(ArgumentError("channel orders differ"))
    check_range(frequencies(a), frequencies(b), "origins")
    epoch_spacing = inv(max(2 * (nsamples(a) - 1), 1) * frequency_spacing(a))
    check_coordinate(epoch(a), epoch(b), epoch_spacing, "epochs", EPOCH_HINT)
    return nothing
end
check_same_grid(::AbstractSeries, ::AbstractSeries) = throw(ArgumentError("domain kinds differ"))

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
    value = metadata_value(boundary)
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
    crop(series::TimeSeries; start, stop)

Copy the half-open interval `[start, stop)` with grid-aligned boundaries. Bounds
are in seconds. Throws
`ArgumentError` for non-finite, off-grid, or out-of-extent boundaries and when
`start` follows `stop`.
"""
function crop(x::TimeSeries; start::Real, stop::Real)
    lo, hi = boundary_index(x, start), boundary_index(x, stop)
    lo <= hi || throw(ArgumentError("crop start must precede stop"))
    return x[(lo + 1):hi, :]
end

"""Zero-pad time sample counts, shifting the origin when padding before."""
function pad(x::TimeSeries; before::Integer = 0, after::Integer = 0)
    before >= 0 && after >= 0 || throw(ArgumentError("padding counts must be nonnegative"))
    n = Base.checked_add(Base.checked_add(nsamples(x), Int(before)), Int(after))
    origin = first(coordinates(x)) - before * step(coordinates(x))
    isfinite(origin) && isfinite(origin + max(n - 1, 0) * step(coordinates(x))) || throw(ArgumentError("padded coordinates must be finite"))
    data = similar(parent(x), (n, size(x, 2)))
    fill!(data, zero(eltype(data)))
    data[(before + 1):(before + nsamples(x)), :] .= parent(x)
    m = merge(metadata(x), (; epoch = origin))
    return rebuild(x, data, m)
end

"""Indices of existing frequency bins inside the inclusive bounds, clipped to the grid."""
function frequency_bins(s::FrequencySeries; fmin = first(frequencies(s)), fmax = max(zero(frequency_spacing(s)), last(frequencies(s))))
    isfinite(fmin) && isfinite(fmax) && fmin <= fmax || throw(ArgumentError("frequency bounds must be finite and ordered"))
    f = frequencies(s)
    return searchsortedfirst(f, fmin):searchsortedlast(f, fmax)
end
