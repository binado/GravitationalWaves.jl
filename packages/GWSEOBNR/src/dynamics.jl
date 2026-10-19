# SPDX-License-Identifier: GPL-3.0-or-later
# Initial-condition and EOM equations from pySEOBNR 0.3.7; see NOTICE.md.
const MTSUN_SI = 4.925490947641267e-6
const MRSUN_SI = 1476.6250380501247
const MPC_SI = 3.085677581491367e22

"""Reject failed solves and evolution ending without a physical termination guard."""
function check_evolution_retcode(retcode::SciMLBase.ReturnCode.T)
    SciMLBase.successful_retcode(retcode) || throw(SolveError("EOB evolution failed: $retcode"))
    retcode == SciMLBase.ReturnCode.Terminated || throw(SolveError("EOB evolution exhausted time coverage without a physical guard"))
    return nothing
end

function check_isco_interval(r_initial, r_final, r_isco)
    r_initial >= r_isco >= r_final || throw(SolveError("Kerr ISCO is unavailable within the evolved trajectory"))
    return nothing
end

function attachment_time(t_isco, delta_t, t_start, t_end)
    t = t_isco - delta_t
    t_start < primal(t) < t_end || throw(SolveError("attachment requires an upstream clamping fallback"))
    return t
end

function eob_rhs(u, pars, t)
    c, coeffs = pars.c, pars.coeffs
    r, phi, pr, L = u
    q, p = (r, phi), (pr, L)
    d = hamiltonian_dynamics(q, p, c)
    omega_circ = orbital_frequency(q, (zero(pr), L), c)
    Fr, Fphi = radiation_reaction(pr, L, d[4], omega_circ, d[5], c, coeffs)
    return SVector(d[6] * d[3], d[4], -d[1] * d[6] + Fr, -d[2] + Fphi)
end

function evolve(m1, m2, chi1, chi2, g::SEOBNRv5HM)
    c = hamiltonian_coefficients(m1, m2, chi1, chi2)
    coeffs = flux_coefficients(c)
    omega0 = π * (m1 + m2) * MTSUN_SI * 0.8g.minimum_frequency
    u0 = initial_conditions(m1, m2, chi1, chi2, omega0)
    af = bbh_final_spin_non_precessing_HBR2016(c.X_1, c.X_2, chi1, chi2, "M3J4")
    r_isco = calc_isco_radius(af)
    ap, am = c.X_1 * chi1 + c.X_2 * chi2, c.X_1 * chi1 - c.X_2 * chi2
    dt_nr = NR_deltaT_NS(c.nu) + NR_deltaT(c.nu, ap, am)
    r_stop = dt_nr > 0 ? 0.98r_isco : 1.4 + zero(r_isco)
    pars = (; c, coeffs)
    condition = function (u, t, integrator)
        primal(u[1]) <= 6 || return false
        d = eob_rhs(u, pars, t)
        omega_rate = dot(hamiltonian_hessian((u[1], u[2]), (u[3], u[4]), c)[4, :], d)
        omega_circ = orbital_frequency((u[1], u[2]), (zero(u[3]), u[4]), c)
        return primal(u[1]) <= primal(r_stop) || primal(d[1]) > 0 || primal(d[3]) > 0 ||
            primal(omega_rate) <= 0 || (primal(u[1]) < 3 && primal(omega_circ) > 1)
    end
    callback = DiscreteCallback(condition, terminate!; save_positions = (true, true))
    # A conservative finite upper bound in geometric units; ending without a guard is failure.
    tmax = 2 * (5 / (256 * primal(c.nu))) * primal(omega0)^(-8 / 3)
    sol = solve(
        ODEProblem(eob_rhs, u0, (0.0, tmax), pars), Vern9();
        reltol = g.reltol, abstol = g.abstol, callback, dense = true,
        dtmax = 10.0, maxiters = 10^7
    )
    check_evolution_retcode(sol.retcode)
    check_isco_interval(sol.u[1][1], sol.u[end][1], r_isco)
    Fisco = t -> sol(t)[1] - r_isco
    t_isco = continuous_root(Fisco, sol.t[1], sol.t[end], t -> sol(t, Val{1})[1])
    t_attach = attachment_time(t_isco, dt_nr, sol.t[1], sol.t[end])
    omega_ref = π * (m1 + m2) * MTSUN_SI * g.reference_frequency
    Fref = t -> orbital_frequency((sol(t)[1], sol(t)[2]), (sol(t)[3], sol(t)[4]), c) - omega_ref
    tref = continuous_root(Fref, sol.t[1], primal(t_attach), t -> dot(hamiltonian_hessian((sol(t)[1], sol(t)[2]), (sol(t)[3], sol(t)[4]), c)[4, :], sol(t, Val{1})))
    return (; sol, c, coeffs, t_isco, t_attach, t_ref = tref, phi_ref = sol(tref)[2], af, r_isco, dt_nr)
end
