# SPDX-License-Identifier: GPL-3.0-or-later
# The original pySEOBNR IC equations are solved as one coupled system;
# see NOTICE.md and test/legacy_initial.jl for the sequential reference.
function dissipative_residual(pr, r, L, c, coeffs)
    q, p = (r, zero(r)), (pr, L)
    hess = hamiltonian_hessian(q, p, c)
    dLdr = -hess[1, 1] / hess[4, 1]
    d = hamiltonian_dynamics(q, p, c)
    omega_circ = orbital_frequency(q, (zero(pr), L), c)
    force = radiation_reaction(pr, L, d[4], omega_circ, d[5], c, coeffs)
    return force[2] / (d[6] * dLdr) - d[3]
end


function initial_residual(u, p)
    r, L, pr = u
    c = hamiltonian_coefficients(p[1], p[2], p[3], p[4])
    d = hamiltonian_gradient((r, zero(r)), (zero(r), L), c)
    dissipative = dissipative_residual(pr, r, L, c, flux_coefficients(c))
    return [d[4] - p[5], d[1], dissipative]
end

function initial_solution(theta; maxiters = 100)
    r = primal(theta[5])^(-2 / 3)
    problem = NonlinearSolve.NonlinearProblem{false}(initial_residual, [r, sqrt(r), -0.001], theta)
    return solve(
        problem, NonlinearSolve.NewtonRaphson(; autodiff = AutoForwardDiff());
        abstol = 1.0e-14, reltol = 1.0e-14, maxiters
    )
end

function check_initial_solution(sol)
    SciMLBase.successful_retcode(sol) || throw(SolveError("NonlinearSolve initial conditions failed: $(sol.retcode)"))
    norm(primal.(sol.resid), Inf) < 1.0e-12 || throw(SolveError("NonlinearSolve initial residual exceeds tolerance"))
    p = primal.(sol.prob.p)
    J = ForwardDiff.jacobian(u -> initial_residual(u, p), primal.(sol.u))
    condition = cond(J)
    isfinite(condition) && condition < 1.0e14 || throw(SolveError("NonlinearSolve initial Jacobian is singular"))
    r, L, pr = sol.u
    primal(r) > 0 && primal(L) > 0 && -0.03 < primal(pr) < 0 || throw(SolveError("NonlinearSolve root is outside the physical initial branch"))
    return nothing
end

function initial_conditions(m1, m2, chi1, chi2, omega)
    theta = collect(promote(m1, m2, chi1, chi2, omega))
    sol = initial_solution(theta)
    check_initial_solution(sol)
    r, L, pr = sol.u
    return SVector(r, zero(r), pr, L)
end
