# SPDX-License-Identifier: GPL-3.0-or-later
# pySEOBNR 0.3.7 (2,2) mode, NQC equations, and merger-ringdown ansatz.
# Modified: continuous evaluation on the ODE interpolant and continuous timing.
function inspiral_mode(e, t)
    u = e.sol(t)
    c, coeffs = e.c, e.coeffs
    q, p = (u[1], u[2]), (u[3], u[4])
    d = hamiltonian_dynamics(q, p, c)
    omega = d[4]
    H = c.nu * d[5]
    omega_circ = orbital_frequency(q, (zero(u[3]), u[4]), c)
    v, vh = omega^(1 / 3), (H * omega)^(1 / 3)
    vs, vhs = ntuple(j -> v^(j - 1), 11), ntuple(j -> vh^(j - 1), 11)
    v_phi2 = (omega / omega_circ^(2 / 3))^2
    newton = coeffs.prefixes[3, 3] * v_phi2 * YLMS_SIGNED[3, 3] * cis(-2 * (u[2] - e.phi_ref))
    source = (H * H - 1.0) / (2.0 * c.nu) + 1.0
    phase = zero(v)
    for j in 1:10
        phase += coeffs.delta[3, 3, j + 1] * vs[j + 1] + coeffs.deltavh[3, 3, j + 1] * vhs[j + 1]
    end
    return complex_tail(2, 2, omega, H) * cis(phase) * source * residual_amplitude(vs, vh, 2, 2, c, coeffs) * newton
end

# Symmetric five-point derivatives of continuous functions in geometric time.
# The spacing is fixed configuration, independent of the parameter draw.
function time_jet(F, t; h = 0.05)
    fm2, fm1, f0, fp1, fp2 = F(t - 2h), F(t - h), F(t), F(t + h), F(t + 2h)
    d1 = (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12h)
    d2 = (-fp2 + 16 * fp1 - 30 * f0 + 16 * fm1 - fm2) / (12h^2)
    return SVector(f0, d1, d2)
end

function nqc_basis(e, t)
    u = e.sol(t)
    omega = orbital_frequency((u[1], u[2]), (u[3], u[4]), e.c)
    r, pr = u[1], u[3]
    p1 = pr / (r * omega)
    q1 = p1^2
    return SVector(q1, q1 / r, q1 / (r * sqrt(r)), p1, p1 * pr^2)
end

function nqc_coefficients(e)
    t = e.t_attach
    0.1 < primal(t) < e.sol.t[end] - 0.1 || throw(SolveError("insufficient NQC matching coverage"))
    amp = s -> abs(inspiral_mode(e, s))
    ampjet = time_jet(amp, t)
    Q = hcat((time_jet(s -> amp(s) * nqc_basis(e, s)[j], t) for j in 1:3)...)
    c = e.c
    target = c.nu * SVector(iv_habs(c), iv_hdot(c), iv_hdotdot(c))
    a = Q \ (target - ampjet)
    # Local phase differences avoid wrapping; the global phase is never unwrapped by samples.
    h0 = inspiral_mode(e, t)
    phasejet = time_jet(s -> angle(inspiral_mode(e, s) / h0), t)
    omega = abs(phasejet[2])
    domega = phasejet[2] * phasejet[3] > 0 ? abs(phasejet[3]) : -abs(phasejet[3])
    P = -hcat((time_jet(s -> nqc_basis(e, s)[j], t)[2:3] for j in 4:5)...)
    b = P \ (SVector(abs(iv_omega(c)), abs(iv_omegadot(c))) - SVector(omega, domega))
    all(isfinite, a) && all(isfinite, b) || throw(SolveError("NQC matching system failed"))
    return (; a, b)
end

function corrected_mode(e, nqc, t)
    basis = nqc_basis(e, t)
    return inspiral_mode(e, t) * (1 + dot(nqc.a, basis[1:3])) * cis(dot(nqc.b, basis[4:5]))
end

function merger_coefficients(e, nqc)
    c, t = e.c, e.t_attach
    Mf = bbh_final_mass_non_precessing_UIB2016(c.X_1, c.X_2, c.chi_1, c.chi_2)
    omega = qnm_220(e.af, Mf)
    sigmaR, sigmaI = imag(omega), -real(omega)
    ampjet = time_jet(s -> abs(corrected_mode(e, nqc, s)), t)
    h0 = corrected_mode(e, nqc, t)
    omega_match = time_jet(s -> angle(corrected_mode(e, nqc, s) / h0), t)[2]
    c1f, c2f = mr_c1f(c), mr_c2f(c)
    if sigmaR > 2c1f * tanh(c2f)
        c1f = sigmaR / (2tanh(c2f))
    end
    d1f, d2f = mr_d1f(c), mr_d2f(c)
    amp, damp = ampjet[1], ampjet[2]
    c1c = (damp - sigmaR * amp) * cosh(c2f)^2 / (c1f * c.nu)
    c2c = amp / c.nu - (damp - sigmaR * amp) * cosh(c2f) * sinh(c2f) / (c1f * c.nu)
    d1c = (omega_match - sigmaI) * (1 + d2f) / (d1f * d2f)
    return (; omega, sigmaR, sigmaI, c1f, c2f, d1f, d2f, c1c, c2c, d1c, h0, eta = c.nu)
end

function ringdown_mode(mr, tau)
    amp = mr.c1c * tanh(mr.c1f * tau + mr.c2f) + mr.c2c
    phase = -mr.d1c * log((1 + mr.d2f * exp(-mr.d1f * tau)) / (1 + mr.d2f))
    return mr.eta * amp * (mr.h0 / abs(mr.h0)) * cis(phase - real(mr.omega) * tau) * exp(mr.sigmaR * tau)
end

function complete_mode(e, nqc, mr, t)
    return t <= e.t_attach ? corrected_mode(e, nqc, t) : ringdown_mode(mr, t - e.t_attach)
end

function peak_time(e, nqc, mr)
    # The (2,2) NQC target has zero amplitude slope at attachment, and the
    # constrained MR ansatz preserves that stationary point. Treat the join
    # explicitly; differentiating across the two pieces would move the peak.
    candidates = [e.t_attach]
    F = t -> time_jet(s -> abs(corrected_mode(e, nqc, s)), t)[2]
    lo = max(e.sol.t[1] + 0.1, primal(e.t_attach) - 100)
    hi = primal(e.t_attach) - 0.15
    grid = range(lo, hi; length = 101)
    for (a, b) in zip(grid[1:(end - 1)], grid[2:end])
        if primal(F(a)) > 0 && primal(F(b)) <= 0
            t = continuous_root(F, a, b, s -> time_jet(x -> abs(corrected_mode(e, nqc, x)), s)[3])
            push!(candidates, t)
        end
    end
    # Analytic MR amplitude slope; its root is solved within its smooth piece.
    R = function (tau)
        z = mr.c1f * tau + mr.c2f
        A = mr.c1c * tanh(z) + mr.c2c
        return mr.c1c * mr.c1f / cosh(z)^2 + mr.sigmaR * A
    end
    if primal(R(0.01)) > 0 && primal(R(100.0)) < 0
        tau = continuous_root(
            R, 0.01, 100.0, t ->
            -2mr.c1c * mr.c1f^2 * tanh(mr.c1f * t + mr.c2f) / cosh(mr.c1f * t + mr.c2f)^2 +
                mr.sigmaR * mr.c1c * mr.c1f / cosh(mr.c1f * t + mr.c2f)^2
        )
        push!(candidates, e.t_attach + tau)
    end
    amplitudes = [primal(abs(complete_mode(e, nqc, mr, t))) for t in candidates]
    return candidates[argmax(amplitudes)]
end
