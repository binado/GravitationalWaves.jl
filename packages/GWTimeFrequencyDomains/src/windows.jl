# Rising Planck transition, with endpoint guards and an overflow-safe logistic.
function planck_transition(x)
    x <= 0 && return zero(x)
    x >= 1 && return one(x)
    z = inv(x) - inv(one(x) - x)
    if z > 0
        e = exp(-z)
        return e / (one(e) + e)
    end
    return inv(one(z) + exp(z))
end

"""
    PlanckWindow(td::TimeDomain; start_width::Real, end_width::Real)

Prepare a Planck taper on the fixed sampling grid `td`. Widths are finite,
nonnegative transition durations in seconds; zero disables the corresponding edge.
Their sum must not exceed the span between stored endpoints, `(nsamples(td)-1)*dt`.
Empty and singleton domains accept only zero widths.

Transitions are anchored at the first and last stored samples and evaluated
continuously, without rounding widths to sample counts. Positive widths shorter
than `dt` are valid. Transitions may meet; a sampled unit plateau need not exist.
Epoch and channels do not affect the weights. The floating element type promotes
spacing and widths to at least `Float64`, preserving `BigFloat` inputs.

The object owns a real `weights` vector and the sample spacing, independent of
epoch and channels. Treat its weights as read-only configuration. Construct it
on the host, outside differentiation, and use `check_compatible(td, window)` and
`check_data(td, h)` there. Apply it using [`apply_window`](@ref) or
[`apply_window!`](@ref) inside numerical kernels.
"""
struct PlanckWindow{D <: Real, T <: AbstractFloat}
    spacing::D
    weights::Vector{T}
    PlanckWindow(spacing, weights, ::Validated) = new{typeof(spacing), eltype(weights)}(spacing, weights)
end

function PlanckWindow(td::TimeDomain; start_width::Real, end_width::Real)
    dt, start, stop = promote(sample_spacing(td), metadata(start_width), metadata(end_width))
    isfinite(start) && start >= 0 || throw(ArgumentError("start width must be finite and nonnegative"))
    isfinite(stop) && stop >= 0 || throw(ArgumentError("end width must be finite and nonnegative"))
    n = nsamples(td)
    span = max(n - 1, 0) * dt
    start <= span && stop <= span - start || throw(ArgumentError("taper widths must fit between stored endpoints"))
    w = Vector{typeof(dt)}(undef, n)
    for i in eachindex(w)
        rising = iszero(start) ? one(dt) : planck_transition(((i - 1) * dt) / start)
        falling = iszero(stop) ? one(dt) : planck_transition(((n - i) * dt) / stop)
        w[i] = rising * falling
    end
    return PlanckWindow(sample_spacing(td), w, Validated())
end

nsamples(w::PlanckWindow) = length(w.weights)
sample_spacing(w::PlanckWindow) = w.spacing

"""Check a time domain's sample count and spacing against a prepared window on the host; ignore epoch and channels."""
function check_compatible(td::TimeDomain, w::PlanckWindow)
    nsamples(td) == nsamples(w) || throw(DimensionMismatch("sample count does not match window"))
    dt, dw = sample_spacing(td), sample_spacing(w)
    check_coordinate(dt, dw, min(dt, dw), "window sample spacings")
    check_coordinate(nsamples(td) * dt, nsamples(w) * dw, min(dt, dw), "window grid extents")
    return nothing
end

"""
    apply_window(window::PlanckWindow, h)

Allocate tapered samples, multiplying each row of a numeric vector or
`(sample, channel)` matrix by the prepared weight. Inputs are unchanged and
ordinary element-type promotion applies. Views and complex samples are supported.
Check data and window compatibility on the host before evaluation; this kernel
does not repeat domain checks. Differentiation is with respect to `h` only.
"""
apply_window(w::PlanckWindow, h::Samples) = w.weights .* h

"""
    apply_window!(h, window::PlanckWindow)
    apply_window!(out, window::PlanckWindow, h)

Apply prepared weights in place, or write tapered `h` into `out`, returning the
destination. Arrays must have matching sample counts and shapes, checked on the
host. Destination storage determines the element type and must represent the
result. The two-argument form modifies `h`; the three-argument form also permits
`out === h`. Keep window weights as read-only configuration.
"""
function apply_window!(out::Samples, w::PlanckWindow, h::Samples)
    out .= w.weights .* h
    return out
end
apply_window!(h::Samples, w::PlanckWindow) = apply_window!(h, w, h)
