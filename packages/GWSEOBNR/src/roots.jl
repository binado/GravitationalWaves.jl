"""Failure of an initial-condition or continuous timing solve."""
struct SolveError <: Exception
    message::String
end
Base.showerror(io::IO, e::SolveError) = print(io, e.message)

function bracket_root(F, lo, hi; xtol = 1.0e-11, maxiters = 100)
    a, b = Float64(primal(lo)), Float64(primal(hi))
    fa, fb = primal(F(a)), primal(F(b))
    all(isfinite, (fa, fb)) || throw(SolveError("nonfinite bracket residual"))
    iszero(fa) && return a
    iszero(fb) && return b
    signbit(fa) != signbit(fb) || throw(SolveError("root is not bracketed"))
    for _ in 1:maxiters
        mid = (a + b) / 2
        fm = primal(F(mid))
        isfinite(fm) || throw(SolveError("nonfinite root residual"))
        (iszero(fm) || b - a <= xtol) && return mid
        if signbit(fm) == signbit(fa)
            a, fa = mid, fm
        else
            b, fb = mid, fm
        end
    end
    throw(SolveError("bracketed root did not converge"))
end

attach_time_partials(t::Real, residual::AbstractFloat, slope) = t
function attach_time_partials(t, residual::ForwardDiff.Dual{Tag, V, P}, slope) where {Tag, V, P}
    s = primal(slope)
    isfinite(s) && abs(s) > 1.0e-14 || throw(SolveError("singular continuous timing root"))
    return ForwardDiff.Dual{Tag}(t, ntuple(j -> -ForwardDiff.partials(residual)[j] / s, P))
end

function continuous_root(F, lo, hi, derivative)
    t = bracket_root(F, lo, hi)
    slope = derivative(t)
    isfinite(primal(slope)) && abs(primal(slope)) > 1.0e-14 || throw(SolveError("singular continuous timing root"))
    return attach_time_partials(t, F(t), slope)
end
