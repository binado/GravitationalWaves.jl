"""Failure of an initial-condition or continuous timing solve."""
struct SolveError <: Exception
    message::String
end
Base.showerror(io::IO, e::SolveError) = print(io, e.message)

function newton_root(F, guess, theta; tolerance = 1.0e-12, maxiters = 60)
    x = collect(guess)
    for _ in 1:maxiters
        residual = F(x, theta)
        all(isfinite, residual) || throw(SolveError("nonfinite initial-condition residual"))
        J = ForwardDiff.jacobian(y -> F(y, theta), x)
        isfinite(cond(J)) && cond(J) < 1.0e14 || throw(SolveError("singular initial-condition Jacobian"))
        norm(residual, Inf) <= tolerance && return x
        dx = J \ residual
        accepted = false
        for j in 0:30
            trial = x - 2.0^(-j) * dx
            trial[1] > 0 || continue
            next_residual = F(trial, theta)
            if all(isfinite, next_residual) && norm(next_residual) < norm(residual)
                x = trial
                accepted = true
                break
            end
        end
        accepted || throw(SolveError("initial-condition Newton line search failed"))
    end
    throw(SolveError("initial-condition root did not converge"))
end

function implicit_root(F, guess, theta::AbstractVector{<:AbstractFloat})
    return newton_root(F, guess, theta)
end

function implicit_root(F, guess, theta::AbstractVector{ForwardDiff.Dual{Tag, V, P}}) where {Tag, V <: AbstractFloat, P}
    t0 = ForwardDiff.value.(theta)
    x = newton_root(F, guess, t0)
    Jx = ForwardDiff.jacobian(y -> F(y, t0), x)
    Jtheta = ForwardDiff.jacobian(y -> F(x, y), t0)
    isfinite(cond(Jx)) && cond(Jx) < 1.0e14 || throw(SolveError("singular implicit root Jacobian"))
    dx = -(Jx \ Jtheta)
    return map(eachindex(x)) do i
        ForwardDiff.Dual{Tag}(
            x[i], ntuple(P) do j
                sum(dx[i, k] * ForwardDiff.partials(theta[k])[j] for k in eachindex(theta))
            end
        )
    end
end

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
