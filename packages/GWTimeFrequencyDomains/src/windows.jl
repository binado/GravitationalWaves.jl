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
    planck_window(td::TimeDomain; start_width::Real, end_width::Real)

Construct a fresh real vector of Planck taper weights on `td`. Widths are finite,
nonnegative transition durations in seconds; zero disables the corresponding edge.
Their sum must not exceed the span between stored endpoints, `(nsamples(td)-1)*dt`.
Empty and singleton domains accept only zero widths.

Transitions are anchored at the first and last stored samples and evaluated
continuously, without rounding widths to sample counts. Positive widths shorter
than `dt` are valid. Transitions may meet; a sampled unit plateau need not exist.
Epoch and channels do not affect the weights. The floating element type promotes
spacing and widths to at least `Float64`, preserving `BigFloat` inputs.

Construct the window on the host, outside differentiation, then apply it with
`w .* h` to a vector or a `(sample, channel)` matrix. Callers are responsible for
matching the window and data grids.
"""
function planck_window(td::TimeDomain; start_width::Real, end_width::Real)
    dt, start, stop = promote(sample_spacing(td), metadata(start_width), metadata(end_width))
    isfinite(start) && start >= 0 || throw(ArgumentError("start width must be finite and nonnegative"))
    isfinite(stop) && stop >= 0 || throw(ArgumentError("end width must be finite and nonnegative"))
    n = nsamples(td)
    span = max(n - 1, 0) * dt
    start <= span && stop <= span - start || throw(ArgumentError("taper widths must fit between stored endpoints"))
    w = Vector{typeof(dt)}(undef, n)
    for i in 0:(n - 1)
        rising = iszero(start) ? one(dt) : planck_transition((i * dt) / start)
        falling = iszero(stop) ? one(dt) : planck_transition(((n - 1 - i) * dt) / stop)
        w[i + 1] = rising * falling
    end
    return w
end
