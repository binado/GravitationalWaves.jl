"""
    NaturalCubicSpline(t, u)

Natural cubic spline (zero second derivative at both ends) through the knots `t` (strictly
increasing reals) with values `u` (real or complex). Mirrors `gwsurrogate`'s
`_spline_interp.cpp`. The second derivatives are solved with the Thomas algorithm in buffers
whose element type is promoted from the inputs, so the construction differentiates with respect
to both `t` and `u`.
"""
struct NaturalCubicSpline{T <: AbstractVector, U <: AbstractVector, M <: AbstractVector}
    t::T
    u::U
    m::M
end

function NaturalCubicSpline(t::AbstractVector{<:Real}, u::AbstractVector{<:Number})
    length(t) == length(u) || throw(DimensionMismatch("t and u must have the same length"))
    knots = SplineKnots(t)
    TU = promote_type(eltype(u), eltype(knots.cp))
    m = zeros(TU, length(t))
    spline_second_derivatives!(m, zeros(TU, length(t)), knots, u)
    return NaturalCubicSpline(t, u, m)
end

# The part of the spline construction that depends on the knots alone: the step sizes and the
# forward-sweep coefficients of the Thomas algorithm for the interior equations
#   h[i-1] m[i-1] + 2 (h[i-1] + h[i]) m[i] + h[i] m[i+1] = 6 ((u[i+1]-u[i])/h[i] - (u[i]-u[i-1])/h[i-1]).
# Splines through the same knots share it.
struct SplineKnots{T}
    h::Vector{T}
    cp::Vector{T}
    denom::Vector{T}
end

function SplineKnots(t::AbstractVector{<:Real})
    n = length(t)
    n >= 3 || throw(ArgumentError("a natural cubic spline needs at least 3 knots"))
    issorted(t; lt = <=) || throw(ArgumentError("knots must be strictly increasing"))
    TT = float(eltype(t))
    h = TT[t[i + 1] - t[i] for i in 1:(n - 1)]
    cp = zeros(TT, n)
    denom = zeros(TT, n)
    for i in 2:(n - 1)
        # cp[1] = 0 starts the recursion, so the first row needs no special case.
        denom[i] = 2 * (h[i - 1] + h[i]) - h[i - 1] * cp[i - 1]
        cp[i] = h[i] / denom[i]
    end
    return SplineKnots{TT}(h, cp, denom)
end

# Second derivatives `m` at the knots for the data `u`. `m` and the scratch `dp` must be zero at the
# ends (they are never written there), so both can be reused across data sets.
function spline_second_derivatives!(m, dp, knots::SplineKnots, u)
    (; h, cp, denom) = knots
    n = length(u)
    for i in 2:(n - 1)
        rhs = 6 * ((u[i + 1] - u[i]) / h[i] - (u[i] - u[i - 1]) / h[i - 1])
        dp[i] = (rhs - h[i - 1] * dp[i - 1]) / denom[i]
    end
    m[n - 1] = dp[n - 1]
    for i in (n - 2):-1:2
        m[i] = dp[i] - cp[i] * m[i + 1]
    end
    return m
end

"""
    spline_eval(spline, x)

Evaluate `spline` at `x`. Outside the knot range the cubic of the first or last segment is
extended.
"""
function spline_eval(s::NaturalCubicSpline, x)
    (; t, u, m) = s
    i = clamp(searchsortedlast(t, x), 1, length(t) - 1)
    h = t[i + 1] - t[i]
    a = (t[i + 1] - x) / h
    b = (x - t[i]) / h
    return a * u[i] + b * u[i + 1] + ((a^3 - a) * m[i] + (b^3 - b) * m[i + 1]) * (h^2 / 6)
end

"""
    spline_columns(t_out, t_in, Y) -> Matrix

Natural cubic spline of every column of `Y` (knots `t_in`), evaluated at `t_out`. Each entry equals
`spline_eval(NaturalCubicSpline(t_in, Y[:, j]), x)`, but the knot factorization and the interval and
weights of each evaluation point are computed once and shared by all columns.
"""
function spline_columns(t_out, t_in, Y::AbstractMatrix)
    n = length(t_in)
    size(Y, 1) == n || throw(DimensionMismatch("Y must have one row per knot"))
    knots = SplineKnots(t_in)
    n_out = length(t_out)
    TW = promote_type(eltype(knots.cp), float(eltype(t_out)))
    interval = Vector{Int}(undef, n_out)
    wa = Vector{TW}(undef, n_out)
    wb = Vector{TW}(undef, n_out)
    wm_lo = Vector{TW}(undef, n_out)
    wm_hi = Vector{TW}(undef, n_out)
    wscale = Vector{TW}(undef, n_out)
    for (k, x) in enumerate(t_out)
        i = clamp(searchsortedlast(t_in, x), 1, n - 1)
        h = knots.h[i]
        a = (t_in[i + 1] - x) / h
        b = (x - t_in[i]) / h
        interval[k] = i
        wa[k] = a
        wb[k] = b
        wm_lo[k] = a^3 - a
        wm_hi[k] = b^3 - b
        wscale[k] = h^2 / 6
    end
    TU = promote_type(eltype(Y), eltype(knots.cp))
    m = zeros(TU, n)
    dp = zeros(TU, n)
    out = Matrix{promote_type(eltype(Y), TW)}(undef, n_out, size(Y, 2))
    for j in axes(Y, 2)
        u = view(Y, :, j)
        spline_second_derivatives!(m, dp, knots, u)
        for k in 1:n_out
            i = interval[k]
            out[k, j] = wa[k] * u[i] + wb[k] * u[i + 1] + (wm_lo[k] * m[i] + wm_hi[k] * m[i + 1]) * wscale[k]
        end
    end
    return out
end
