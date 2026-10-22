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
    PlanckWindow(td::TimeSeries; start_width::Real, end_width::Real)

Prepare a Planck taper on the fixed sampling grid `td`. Widths are finite,
nonnegative transition durations in seconds; zero disables the corresponding edge.
Their sum must not exceed the span between stored endpoints, `(nsamples(td)-1)*dt`.
Empty and singleton series accept only zero widths.

Transitions are anchored at the first and last stored samples and evaluated
continuously, without rounding widths to sample counts. Positive widths shorter
than `dt` are valid. Transitions may meet; a sampled unit plateau need not exist.
Epoch and channels do not affect the weights. The floating element type promotes
spacing and widths to at least `Float64`, preserving `BigFloat` inputs.

The object owns a real `weights` vector and the sample spacing, independent of
epoch and channels. Treat its weights as read-only configuration. Construct it
on the host, outside differentiation, and use `check_compatible(td, window)` there.
Use `check_data(td, h)` for vector/matrix layouts. Apply it using [`apply_window`](@ref) or
[`apply_window!`](@ref) inside numerical kernels.
"""
struct PlanckWindow{D <: Real, T <: AbstractFloat}
    spacing::D
    weights::Vector{T}
    PlanckWindow(spacing, weights, ::Validated) = new{typeof(spacing), eltype(weights)}(spacing, weights)
end

function PlanckWindow(td::TimeSeries; start_width::Real, end_width::Real)
    dt, start, stop = promote(step(td), metadata_value(start_width), metadata_value(end_width))
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
    return PlanckWindow(step(td), w, Validated())
end

nsamples(w::PlanckWindow) = length(w.weights)
Base.step(w::PlanckWindow) = w.spacing

"""Check a time series sample count and spacing against a prepared window on the host; ignore epoch and channels."""
function check_compatible(td::TimeSeries, w::PlanckWindow)
    nsamples(td) == nsamples(w) || throw(DimensionMismatch("sample count does not match window"))
    dt, dw = step(td), step(w)
    check_coordinate(dt, dw, min(dt, dw), "window sample spacings")
    check_coordinate(nsamples(td) * dt, nsamples(w) * dw, min(dt, dw), "window grid extents")
    return nothing
end

"""
    apply_window(window::PlanckWindow, h)

Allocate tapered samples by broadcasting the weights with a numeric array `h`.
Weights act along the first dimension; trailing dimensions broadcast normally.
Inputs are unchanged and
ordinary element-type promotion applies. Views and complex samples are supported.
Check window/grid compatibility on the host before evaluation; for vector/matrix
layouts, use `check_data` for shape validation. Higher-dimensional arrays
are also supported when their axes broadcast with the weights. This kernel does
not repeat series checks. Differentiation is with respect to `h` only.
"""
apply_window(w::PlanckWindow, h::AbstractArray{<:Number}) = w.weights .* h

"""
    apply_window!(h, window::PlanckWindow)
    apply_window!(out, window::PlanckWindow, h)

Apply prepared weights in place, or write tapered `h` into `out`, returning the
destination. Weights and `h` must broadcast to the axes of `out`. Destination
storage determines the element type and must represent the
result. The two-argument form modifies `h`; the three-argument form also permits
`out === h`. Keep window weights as read-only configuration.
"""
function apply_window!(out::AbstractArray{<:Number}, w::PlanckWindow, h::AbstractArray{<:Number})
    out .= w.weights .* h
    return out
end
apply_window!(h::AbstractArray{<:Number}, w::PlanckWindow) = apply_window!(h, w, h)

function PlanckWindow(; nsamples::Integer, sampling_interval::Real, start_width::Real, end_width::Real)
    nsamples >= 0 || throw(ArgumentError("sample count must be nonnegative"))
    # A range is a storage-free prototype; no sample buffer is allocated.
    prototype = TimeSeries(range(0.0; length = nsamples, step = 0.0); epoch = 0.0, sampling_interval, channels = (:window,))
    return PlanckWindow(prototype; start_width, end_width)
end
function apply_window(w::PlanckWindow, h::TimeSeries)
    check_compatible(h, w)
    return rebuild(h, apply_window(w, parent(h)), metadata(h))
end
function apply_window!(out::TimeSeries, w::PlanckWindow, h::TimeSeries)
    check_same_grid(out, h)
    check_compatible(h, w)
    apply_window!(parent(out), w, parent(h))
    return out
end
apply_window!(h::TimeSeries, w::PlanckWindow) = apply_window!(h, w, h)
