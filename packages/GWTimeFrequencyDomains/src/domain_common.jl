abstract type AbstractDomain end
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

"""Segment epoch in seconds, independent of sample precision."""
epoch(d::AbstractDomain) = d.epoch
"""Ordered, unique Symbol channel names."""
channels(d::AbstractDomain) = d.channels

"""Check numeric vector/matrix storage, one-based axes, sample count and channel count on the host."""
function check_data(d::AbstractDomain, values::AbstractVecOrMat{<:Number})
    Base.require_one_based_indexing(values)
    size(values, 1) == nsamples(d) || throw(DimensionMismatch("sample count does not match domain"))
    size(values, 2) == length(channels(d)) || throw(DimensionMismatch("channel count does not match domain"))
    return nothing
end
check_data(::AbstractDomain, values) = throw(ArgumentError("data must be a numeric vector or matrix, got $(typeof(values))"))

"""Host compatibility check extended by packages owning prepared numerical objects."""
function check_compatible end

"""Check even-length one-sided real-FFT eligibility on the host."""
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
