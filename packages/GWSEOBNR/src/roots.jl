"""Failure of an initial-condition or continuous timing solve."""
struct SolveError <: Exception
    message::String
end
Base.showerror(io::IO, e::SolveError) = print(io, e.message)

function bracket_root(F, lo, hi, args...; xtol = 1.0e-11, maxiters = 100)
    a, b = Float64(primal(lo)), Float64(primal(hi))
    fa, fb = Float64(primal(stopgrad(F(a, args...)))), Float64(primal(stopgrad(F(b, args...))))
    all(isfinite, (fa, fb)) || throw(SolveError("nonfinite bracket residual"))
    iszero(fa) && return a
    iszero(fb) && return b
    (fa < 0) != (fb < 0) || throw(SolveError("root is not bracketed"))
    for _ in 1:maxiters
        mid = (a + b) / 2
        fm = Float64(primal(stopgrad(F(mid, args...))))
        isfinite(fm) || throw(SolveError("nonfinite root residual"))
        (iszero(fm) || b - a <= xtol) && return mid
        if (fm < 0) == (fa < 0)
            a, fa = mid, fm
        else
            b, fb = mid, fm
        end
    end
    throw(SolveError("bracketed root did not converge"))
end

# Value-preserving stop-gradient: the Enzyme extension marks it inactive.
stopgrad(x::Real) = x
stopgrad(x::ForwardDiff.Dual) = ForwardDiff.value(x)

# Implicit-function correction at a converged root t0 of F(t, theta) = 0:
#   t = t0 - (F(t0) - stopgrad(F(t0))) / F'(t0),
# `F(t, args...)` and `derivative(t, args...)` take their data as arguments
# rather than captured variables: Enzyme reverse mode handles closures poorly.
# so t == t0 while dt/dtheta = -dF/dtheta / F'. Bisection sees only
# stop-gradient float residuals and the slope is a constant first-order factor, so this
# expresses the implicit pullback with ordinary AD primitives for any backend.
function continuous_root(F, lo, hi, derivative, args...)
    t0 = bracket_root(F, lo, hi, args...)
    slope = stopgrad(derivative(t0, args...))
    isfinite(primal(slope)) && abs(primal(slope)) > 1.0e-14 || throw(SolveError("singular continuous timing root"))
    residual = F(t0, args...)
    return t0 - (residual - stopgrad(residual)) / slope
end
