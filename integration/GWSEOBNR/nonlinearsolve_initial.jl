# Experimental initial-condition solver; production waveforms do not call it.
module NonlinearInitialPrototype

using GWSEOBNR, NonlinearSolve, ForwardDiff, Enzyme, SciMLSensitivity
using LinearAlgebra
include("../../packages/GWSEOBNR/test/legacy_initial.jl")
const S = GWSEOBNR
const DEFAULT_SENSEALG = SteadyStateAdjoint(; autojacvec = false)

# Couple the original two conservative equations to the original dissipative
# equation. The conservative equations do not depend on radial momentum, so
# this has the same roots as the production solver's sequential construction.
function initial_residual(u, p)
    r, L, pr = u
    c = S.hamiltonian_coefficients(p[1], p[2], p[3], p[4])
    d = S.hamiltonian_gradient((r, zero(r)), (zero(r), L), c)
    dissipative = S.dissipative_residual(pr, r, L, c, S.flux_coefficients(c))
    return [d[4] - p[5], d[1], dissipative]
end

function initial_solution(x; maxiters = 100, sensealg = DEFAULT_SENSEALG)
    omega = π * (x[1] + x[2]) * S.MTSUN_SI * 16
    p = [x[1], x[2], x[3], x[4], omega]
    r = S.primal(omega)^(-2 / 3)
    problem = NonlinearProblem{false}(initial_residual, [r, sqrt(r), -0.001], p)
    return solve(
        problem, NewtonRaphson(; autodiff = NonlinearSolve.ADTypes.AutoForwardDiff());
        abstol = 1.0e-14, reltol = 1.0e-14, maxiters, sensealg
    )
end

function check_initial_solution(sol)
    S.SciMLBase.successful_retcode(sol) || throw(S.SolveError("NonlinearSolve initial conditions failed: $(sol.retcode)"))
    norm(S.primal.(sol.resid), Inf) < 1.0e-12 || throw(S.SolveError("NonlinearSolve initial residual exceeds tolerance"))
    p = S.primal.(sol.prob.p)
    J = ForwardDiff.jacobian(u -> initial_residual(u, p), S.primal.(sol.u))
    condition = cond(J)
    isfinite(condition) && condition < 1.0e14 || throw(S.SolveError("NonlinearSolve initial Jacobian is singular"))
    r, L, pr = sol.u
    S.primal(r) > 0 && S.primal(L) > 0 && -0.03 < S.primal(pr) < 0 || throw(S.SolveError("NonlinearSolve root is outside the physical initial branch"))
    return nothing
end

# Validation executes in the primal; the solved state remains active.
Enzyme.EnzymeRules.inactive(::typeof(check_initial_solution), sol) = true

function initial_conditions(x; sensealg = DEFAULT_SENSEALG)
    sol = initial_solution(x; sensealg)
    check_initial_solution(sol)
    u = sol.u
    return [u[1], zero(u[1]), u[3], u[2]]
end

function reference_initial(x)
    omega = π * (x[1] + x[2]) * S.MTSUN_SI * 16
    return collect(LegacyInitial.initial_conditions(x..., omega))
end

function initial_loss(x; sensealg = DEFAULT_SENSEALG)
    u = initial_conditions(x; sensealg)
    return u[1] + 3u[3] + 2u[4]
end
reference_loss(x) = (u = reference_initial(x); u[1] + 3u[3] + 2u[4])

function fivepoint(f, x, h)
    return map(eachindex(x)) do j
        v = zeros(length(x)); v[j] = h
        (-f(x + 2v) + 8f(x + v) - 8f(x - v) + f(x - 2v)) / (12h)
    end
end

relative(a, b) = norm(a - b) / max(norm(b), eps(Float64))

end
