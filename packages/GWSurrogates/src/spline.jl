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
    n = length(t)
    n == length(u) || throw(DimensionMismatch("t and u must have the same length"))
    n >= 3 || throw(ArgumentError("a natural cubic spline needs at least 3 knots"))
    issorted(t; lt = <=) || throw(ArgumentError("knots must be strictly increasing"))
    TT = float(eltype(t))
    TU = promote_type(eltype(u), TT)
    m = zeros(TU, n)
    # Thomas algorithm for the interior equations
    #   h[i-1] m[i-1] + 2 (h[i-1] + h[i]) m[i] + h[i] m[i+1] = 6 ((u[i+1]-u[i])/h[i] - (u[i]-u[i-1])/h[i-1]).
    cp = zeros(TT, n)
    dp = zeros(TU, n)
    for i in 2:(n - 1)
        hl = t[i] - t[i - 1]
        hr = t[i + 1] - t[i]
        rhs = 6 * ((u[i + 1] - u[i]) / hr - (u[i] - u[i - 1]) / hl)
        if i == 2
            denom = 2 * (hl + hr)
            cp[i] = hr / denom
            dp[i] = rhs / denom
        else
            denom = 2 * (hl + hr) - hl * cp[i - 1]
            cp[i] = hr / denom
            dp[i] = (rhs - hl * dp[i - 1]) / denom
        end
    end
    m[n - 1] = dp[n - 1]
    for i in (n - 2):-1:2
        m[i] = dp[i] - cp[i] * m[i + 1]
    end
    return NaturalCubicSpline(t, u, m)
end

"""
    spline_eval(spline, x)

Evaluate `spline` at `x`. Outside the knot range the cubic of the first or last segment is
extended.
"""
function spline_eval(s::NaturalCubicSpline, x)
    t = s.t
    i = clamp(searchsortedlast(t, x), 1, length(t) - 1)
    h = t[i + 1] - t[i]
    a = (t[i + 1] - x) / h
    b = (x - t[i]) / h
    return a * s.u[i] + b * s.u[i + 1] + ((a^3 - a) * s.m[i] + (b^3 - b) * s.m[i + 1]) * (h^2 / 6)
end
