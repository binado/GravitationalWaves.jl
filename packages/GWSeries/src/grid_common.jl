abstract type AbstractSeries{T} <: AbstractMatrix{T} end
struct Validated end

"""Promote grid metadata to at least `Float64`; wider types such as `BigFloat` are preserved."""
metadata_value(x::Real) = float(promote(x, 0.0)[1])

"""Validate channel names (nonempty, unique) and return them as a tuple of `Symbol`s."""
function channel_order(names)
    labels = Tuple(names)
    !isempty(labels) && all(x -> x isa Symbol || x isa AbstractString, labels) ||
        throw(ArgumentError("channels must be nonempty symbol or string names"))
    canonical = map(Symbol, labels)
    allunique(canonical) || throw(ArgumentError("channel names must be unique"))
    return canonical
end

"""Check numeric vector/matrix storage, one-based axes, sample count and channel count on the host."""
function check_data(d::AbstractSeries, values::AbstractVecOrMat{<:Number})
    Base.require_one_based_indexing(values)
    size(values, 1) == nsamples(d) || throw(DimensionMismatch("sample count does not match the series"))
    size(values, 2) == length(channels(d)) || throw(DimensionMismatch("channel count does not match the series"))
    return nothing
end
check_data(::AbstractSeries, values) = throw(ArgumentError("data must be a numeric vector or matrix, got $(typeof(values))"))

"""Host compatibility check extended by packages owning prepared numerical objects."""
function check_compatible end

"""Check real-FFT eligibility on the host; frequency input requires explicit nsamples."""
function check_fourier_domain end

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

"""Throw `ArgumentError` unless coordinates `a` and `b` agree within the few-ulp grid tolerance for `spacing`."""
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
