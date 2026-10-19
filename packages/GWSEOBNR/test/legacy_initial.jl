# SPDX-License-Identifier: GPL-3.0-or-later
# Frozen sequential initial-condition solver from commit 6bc5a32.
# Test/reference tooling only; physics helpers remain tied to GWSEOBNR's port.
module LegacyInitial
using GWSEOBNR, ForwardDiff, LinearAlgebra
using StaticArrays: SVector
using GWSEOBNR: primal, hamiltonian_coefficients, hamiltonian_gradient,
    dissipative_residual, flux_coefficients
const SolveError = GWSEOBNR.SolveError

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

function conservative_residual(x, theta)
    c = hamiltonian_coefficients(theta[1], theta[2], theta[3], theta[4])
    d = hamiltonian_gradient((x[1], zero(x[1])), (zero(x[1]), x[2]), c)
    return [d[4] - theta[5], d[1]]
end

function initial_conditions(m1, m2, chi1, chi2, omega)
    theta = collect(promote(m1, m2, chi1, chi2, omega))
    r_guess = primal(omega)^(-2 / 3)
    r, L = implicit_root(conservative_residual, [r_guess, sqrt(r_guess)], theta)
    c = hamiltonian_coefficients(theta[1], theta[2], theta[3], theta[4])
    coeffs = flux_coefficients(c)
    F = pr -> dissipative_residual(pr, r, L, c, coeffs)
    root = bracket_root(F, -0.03, 0.0; xtol = 1.0e-13)
    slope = ForwardDiff.derivative(
        pr -> dissipative_residual(
            pr, primal(r), primal(L),
            hamiltonian_coefficients(primal(m1), primal(m2), primal(chi1), primal(chi2)),
            flux_coefficients(hamiltonian_coefficients(primal(m1), primal(m2), primal(chi1), primal(chi2)))
        ), root
    )
    pr = attach_time_partials(root, F(root), slope)
    return SVector(r, zero(r), pr, L)
end

function initial_conditions(x::AbstractVector)
    omega = π * (x[1] + x[2]) * GWSEOBNR.MTSUN_SI * 16
    return collect(initial_conditions(x..., omega))
end
end
