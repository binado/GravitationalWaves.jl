# IMRPhenomXAS default aligned-spin quadrupole model (Pratten et al. 2020).
# GPL-2.0-or-later, following LALSimulation and the pinned Ripple fit tables.

const XAS_G4 = (0.0, 0.25, 0.75, 1.0)
const XAS_G5 = (0.0, 0.5 - 1 / (2sqrt(2.0)), 0.5, 0.5 + 1 / (2sqrt(2.0)), 1.0)
const XAS_DPHASE0 = 5 / (128 * Float64(π)^(5 / 3))

# The tiny collocation systems retain pivoting and differentiate both A and b.
function xas_solve(rows::NTuple{N, Any}, b::NTuple{N, Any}) where {N}
    A = SMatrix{N, N}(ntuple(k -> rows[mod1(k, N)][cld(k, N)], Val(N * N)))
    return Tuple(lu(A) \ SVector(b))
end

function xas_isco(a)
    r = if abs(a) < 1.0e-4
        # Regular implicit ISCO equation removes sign(a)*sqrt(a^2) at a=0.
        value = 6 - 4sqrt(2 / 3) * a - (7 / 18) * a * a
        for _ in 1:3
            value -= (value * value - 6value + 8a * sqrt(value) - 3a * a) / (2value - 6 + 4a / sqrt(value))
        end
        value
    else
        z1 = 1 + cbrt(1 - a * a) * (cbrt(1 + a) + cbrt(1 - a))
        z1 = min(z1, oftype(z1, 3))
        z2 = sqrt(3a * a + z1 * z1)
        3 + z2 - sign(a) * sqrt((3 - z1) * (3 + z1 + 2z2))
    end
    return inv((sqrt(r)^3 + a) * Float64(π))
end

function xas_phase_cv(i, r, S)
    c = XAS_PHASE_FITS[i]
    return xas_nospin_cv(ntuple(j -> c[j], Val(10)), r.eta) +
        xas_eqspin_cv(ntuple(j -> c[j + 10], Val(29)), r.eta, S) +
        xas_unequal_spin_cv(ntuple(j -> c[j + 39], Val(10)), r.eta, S, r.chia, r.delta)
end

function xas_amp_cv(i, r, S)
    c = XAS_AMP_FITS[i]
    return xas_amp_nospin_cv(ntuple(j -> c[j], Val(8)), r.eta) +
        xas_amp_eqspin_cv(ntuple(j -> c[j + 8], Val(28)), r.eta, S) +
        xas_amp_unequal_spin_cv(ntuple(j -> c[j + 36], Val(6)), r.eta, S, r.chia, r.delta)
end

function xas_inspiral_phase(coefficients::NamedTuple, f)
    pn, logs, sigma = coefficients.pn, coefficients.logs, coefficients.sigma
    phi = sum(ntuple(i -> pn[i] * f^((i - 1) / 3), Val(9)))
    phi += logs[1] * f^(5 / 3) * log(f) + logs[2] * f^2 * log(f) + logs[3] * f^(8 / 3) * log(f)
    phi += sum(ntuple(i -> sigma[i] * f^((i + 7) / 3), Val(4)))
    return -(3 * Float64(π)^(-5 / 3) / 128) * phi * f^(-5 / 3)
end

function xas_dinspiral_phase(coefficients::NamedTuple, f)
    pn, logs, sigma = coefficients.pn, coefficients.logs, coefficients.sigma
    value = sum(ntuple(i -> pn[i] * ((i - 6) / 3) * f^((i - 9) / 3), Val(9)))
    value += logs[1] / f + logs[2] * f^(-2 / 3) * (1 + log(f) / 3) + logs[3] * (1 + log(f))
    value += sum(ntuple(i -> sigma[i] * ((i + 2) / 3) * f^((i - 1) / 3), Val(4)))
    return -(3 * Float64(π)^(-5 / 3) / 128) * value
end

function xas_rd_phase(coefficients::NamedTuple, f)
    c = coefficients.c
    return c[1] * f + 1.5c[2] * f^(2 / 3) - c[3] / f - c[4] / (3f^3) +
        coefficients.c_l / coefficients.fdamp * atan((f - coefficients.f_rd) / coefficients.fdamp)
end

function xas_drd_phase(coefficients::NamedTuple, f)
    c = coefficients.c
    return c[1] + c[2] * f^(-1 / 3) + c[3] / f^2 + c[4] / f^4 + coefficients.c_l / (coefficients.fdamp^2 + (f - coefficients.f_rd)^2)
end

function xas_int_phase(coefficients::NamedTuple, f)
    b = coefficients.b
    return b[1] * f + b[2] * log(f) - b[3] / f - b[4] / (2f^2) - b[5] / (3f^3) +
        2coefficients.c_l / coefficients.fdamp * atan((f - coefficients.f_rd) / (2coefficients.fdamp))
end

function xas_dint_phase(coefficients::NamedTuple, f)
    b = coefficients.b
    return b[1] + b[2] / f + b[3] / f^2 + b[4] / f^3 + b[5] / f^4 + 4coefficients.c_l / (4coefficients.fdamp^2 + (f - coefficients.f_rd)^2)
end

function xas_phase_coefficients(r, chi1, chi2)
    pn = xas_pn_phase(r.eta, r.delta, chi1, chi2)
    cp = ntuple(i -> 0.0026 + XAS_G4[i] * (1.02r.f_meco - 0.0026), Val(4))
    cv = ntuple(i -> xas_phase_cv(i, r, r.S), Val(4))
    rhs = (cv[1] + cv[3], cv[2] + cv[3], cv[3], cv[4] + cv[3])
    insc = xas_solve(ntuple(i -> (one(cp[i]), cp[i]^(1 / 3), cp[i]^(2 / 3), cp[i]), Val(4)), rhs)
    ins = (; pn = pn[1:9], logs = pn[10:12], sigma = ntuple(i -> -5 / (i + 2) * insc[i], Val(4)))
    fm = 0.6 * (0.5r.f_rd + r.f_isco)
    f1 = r.f_meco - 0.03 * (fm - r.f_meco)
    f2 = fm + 0.015 * (fm - r.f_meco)
    rdcp0 = ntuple(i -> fm + XAS_G5[i] * (r.f_rd + 1.25r.fdamp - fm), Val(5))
    rdcp = (rdcp0[1], rdcp0[2], rdcp0[3], r.f_rd, rdcp0[5])
    v = ntuple(i -> xas_phase_cv(i + 8, r, r.s_tot_r), Val(5))
    rdvalues = (v[1] + v[2] + v[4], v[2] + v[4], v[3] + v[4], v[4], v[5] + v[4])
    rdc = xas_solve(
        ntuple(
            i -> (
                one(rdcp[i]), rdcp[i]^(-1 / 3), rdcp[i]^-2,
                rdcp[i]^-4, -XAS_DPHASE0 / (r.fdamp^2 + (rdcp[i] - r.f_rd)^2),
            ), Val(5)
        ), rdvalues
    )
    c_l = -XAS_DPHASE0 * rdc[5]
    rd = (; c = rdc[1:4], c_l, f_rd = r.f_rd, fdamp = r.fdamp)
    vi = ntuple(i -> xas_phase_cv(i + 4, r, r.s_tot_r), Val(4))
    intvalues = (
        xas_dinspiral_phase(ins, f1), 0.75 * (vi[1] + v[4]) + 0.25vi[3],
        vi[2] + v[4], vi[4] + vi[2] + v[4], rdvalues[1],
    )
    icp = ntuple(i -> f1 + XAS_G5[i] * (f2 - f1), Val(5))
    irhs = ntuple(i -> intvalues[i] - 4c_l / (4r.fdamp^2 + (icp[i] - r.f_rd)^2), Val(5))
    bc = xas_solve(ntuple(i -> ntuple(j -> (r.f_rd / icp[i])^(j - 1), Val(5)), Val(5)), irhs)
    int = (; b = ntuple(i -> bc[i] * r.f_rd^(i - 1), Val(5)), c_l, f_rd = r.f_rd, fdamp = r.fdamp)
    alpha1 = xas_dinspiral_phase(ins, f1) - xas_dint_phase(int, f1)
    alpha0 = xas_inspiral_phase(ins, f1) - xas_int_phase(int, f1) - alpha1 * f1
    beta1 = xas_dint_phase(int, f2) + alpha1 - xas_drd_phase(rd, f2)
    beta0 = xas_int_phase(int, f2) + alpha0 + alpha1 * f2 - xas_rd_phase(rd, f2) - beta1 * f2
    return (; ins, int, rd, f1, f2, alpha0, alpha1, beta0, beta1, eta = r.eta)
end

function phenomxas_phase(coefficients::NamedTuple, f)
    value = if f < coefficients.f1
        xas_inspiral_phase(coefficients.ins, f)
    elseif f > coefficients.f2
        xas_rd_phase(coefficients.rd, f) + coefficients.beta0 + coefficients.beta1 * f
    else
        xas_int_phase(coefficients.int, f) + coefficients.alpha0 + coefficients.alpha1 * f
    end
    return value / coefficients.eta
end

function xas_dphase(coefficients::NamedTuple, f)
    value = if f < coefficients.f1
        xas_dinspiral_phase(coefficients.ins, f)
    elseif f > coefficients.f2
        xas_drd_phase(coefficients.rd, f) + coefficients.beta1
    else
        xas_dint_phase(coefficients.int, f) + coefficients.alpha1
    end
    return value / coefficients.eta
end

function xas_inspiral_amp(coefficients::NamedTuple, f)
    return sum(ntuple(i -> coefficients.pn[i] * f^((i - 1) / 3), Val(7))) +
        sum(ntuple(i -> coefficients.rho[i] * f^((i + 6) / 3), Val(3)))
end

function xas_dinspiral_amp(coefficients::NamedTuple, f)
    return sum(ntuple(i -> coefficients.pn[i + 1] * (i / 3) * f^((i - 3) / 3), Val(6))) +
        sum(ntuple(i -> coefficients.rho[i] * ((i + 6) / 3) * f^((i + 3) / 3), Val(3)))
end

xas_rd_amp(coefficients::NamedTuple, f) = exp(-(f - coefficients.f_rd) * coefficients.gamma_r) * coefficients.gamma_d13 / ((f - coefficients.f_rd)^2 + coefficients.gamma_d2)
xas_drd_amp(coefficients::NamedTuple, f) = xas_rd_amp(coefficients, f) * (-coefficients.gamma_r - 2(f - coefficients.f_rd) / ((f - coefficients.f_rd)^2 + coefficients.gamma_d2))

function xas_amplitude_coefficients(r, chi1, chi2)
    f1 = r.f_meco + 0.25 * (r.f_isco - r.f_meco)
    cv = ntuple(i -> xas_amp_cv(i, r, r.S), Val(3))
    cp = (0.5 * f1, 0.75 * f1, f1)
    rho = xas_solve(ntuple(i -> ntuple(j -> cp[i]^((j + 6) / 3), Val(3)), Val(3)), cv)
    ins = (; pn = xas_pn_amplitude(r.eta, r.delta, chi1, chi2), rho)
    gamma2 = xas_amp_cv(5, r, r.s_tot_r)
    gamma3 = xas_amp_cv(6, r, r.s_tot_r)
    shift = gamma2 <= 1 ? (sqrt(1 - gamma2^2) - 1) / gamma2 : -inv(gamma2)
    f4 = abs(r.f_rd + r.fdamp * gamma3 * shift)
    vrd = xas_amp_cv(7, r, r.s_tot_r)
    gamma1 = vrd / (r.fdamp * gamma3) * ((f4 - r.f_rd)^2 + (r.fdamp * gamma3)^2) * exp((f4 - r.f_rd) * gamma2 / (r.fdamp * gamma3))
    rd = (; f_rd = r.f_rd, gamma_r = gamma2 / (r.fdamp * gamma3), gamma_d2 = (gamma3 * r.fdamp)^2, gamma_d13 = r.fdamp * gamma1 * gamma3)
    a1, a4 = xas_inspiral_amp(ins, f1), xas_rd_amp(rd, f4)
    v1, v4 = f1^(7 / 6) / a1, f4^(7 / 6) / a4
    v2 = inv(xas_amp_cv(4, r, r.s_tot_r))
    d1 = (7 / 6) * f1^(1 / 6) / a1 - f1^(7 / 6) * xas_dinspiral_amp(ins, f1) / a1^2
    d4 = (7 / 6) * f4^(1 / 6) / a4 - f4^(7 / 6) * xas_drd_amp(rd, f4) / a4^2
    width = f4 - f1
    # Reciprocal amplitude quartic in t=(f-f1)/width, with endpoint slopes.
    delta = xas_solve(
        (
            (1.0, 0.0, 0.0, 0.0, 0.0), (1.0, 0.5, 0.25, 0.125, 0.0625),
            (1.0, 1.0, 1.0, 1.0, 1.0), (0.0, 1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 2.0, 3.0, 4.0),
        ),
        (v1, v2, v4, d1 * width, d4 * width)
    )
    return (; ins, rd, f1, f4, width, delta)
end

function phenomxas_amplitude(a, f)
    value = if f < a.f1
        xas_inspiral_amp(a.ins, f)
    elseif f > a.f4
        xas_rd_amp(a.rd, f)
    else
        t = (f - a.f1) / a.width
        f^(7 / 6) / evalpoly(t, a.delta)
    end
    return value * f^(-7 / 6)
end

function phenomxas_model(m1in, m2in, chi1in, chi2in, fref, phase)
    m1, m2, chi1, chi2 = m1in, m2in, chi1in, chi2in
    r = xas_remnant(m1, m2, chi1, chi2)
    amplitude = xas_amplitude_coefficients(r, chi1, chi2)
    phase_coefficients = xas_phase_coefficients(r, chi1, chi2)
    M_sec = (m1 + m2) * MTSUN_SI
    _, linb, psi4 = xas_phase_at_peak(r.eta, r.s_tot_r, r.chia, r.delta)
    linb -= xas_dphase(phase_coefficients, r.f_rd - r.fdamp) + 2Float64(π) * (500 + psi4)
    Mfref = M_sec * fref
    # ChooseFDWaveform's polarization convention includes the additional π.
    phiref = -phenomxas_phase(phase_coefficients, Mfref) - linb * Mfref + 2phase + PI_4 - Float64(π)
    return (; r, amplitude, phase = phase_coefficients, M_sec, linb, phiref)
end

# IMRPhenomXAS on the supplied frequencies: nonnegative amplitude and unwrapped phase, zero
# outside `[f_start, f_max]` and above the inclusive cutoff `0.3 / M_sec`. The inclination
# factors are applied by `hphc_fd`. Inputs are trusted; run `validate_parameters` on the host first.
function amplitude_phase_fd(f::AbstractRange{Float64}, params::NamedTuple, g::IMRPhenomXAS)
    T, p = promote_parameters(aligned_parameters(params))
    f_min, f_max = frequency_bounds(g, f)
    f_start = start_frequency(g, f_min)
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    lo, hi = phenom_frequency_slice(f, f_start, f_min, f_max)
    lo > hi && return (; amplitude, phase)

    f_ref = phenom_reference_frequency(g, f_start)
    model = phenomxas_model(p.mass_1, p.mass_2, p.chi_1, p.chi_2, f_ref, p.phase)
    amp0 = phenom_amplitude_scale(p.mass_1, p.mass_2, p.luminosity_distance)
    norm = amp0 * sqrt(2model.r.eta / 3) * Float64(π)^(-1 / 6)
    cutoff = 0.3 / model.M_sec
    for i in lo:hi
        amplitude[i], phase[i] = phenomxas_amplitude_phase(f[i], model, norm, cutoff)
    end
    return (; amplitude, phase)
end

# Scalar evaluation at a positive frequency in Hz, using coefficients prepared once per
# waveform; the caller handles the lower bounds.
function phenomxas_amplitude_phase(frequency, model, norm, cutoff)
    frequency <= cutoff || return (zero(norm), zero(norm))
    Mf = model.M_sec * frequency
    amp = norm * phenomxas_amplitude(model.amplitude, Mf)
    phi = phenomxas_phase(model.phase, Mf) + model.linb * Mf + model.phiref
    return canonical_amplitude_phase(amp, phi)
end

function check_waveform_inputs(g::IMRPhenomXAS, delta_f, params::NamedTuple, f_start, f_max)
    (delta_f > 0 && f_start > 0 && f_max > f_start && g.reference_frequency >= 0) || throw(ArgumentError("invalid frequency bounds for IMRPhenomXAS"))
    check_phenom_parameters(params)
    params.mass_1 >= params.mass_2 || throw(ArgumentError("IMRPhenomXAS requires mass_1 ≥ mass_2"))
    params.mass_1 / params.mass_2 <= 1000 || throw(ArgumentError("IMRPhenomXAS requires mass ratio <= 1000"))
    r = xas_remnant(params.mass_1, params.mass_2, params.chi_1, params.chi_2)
    (abs(r.a) < 1 && 0 < r.Erad < 1 && r.f_meco <= r.f_isco) || throw(ArgumentError("invalid remnant or MECO/ISCO frequencies"))
    0.3 / ((params.mass_1 + params.mass_2) * MTSUN_SI) > f_start || throw(ArgumentError("XAS cutoff must exceed the starting frequency"))
    return nothing
end
