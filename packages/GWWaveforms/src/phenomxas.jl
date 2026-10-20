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
    return xas_nospin_CV(ntuple(j -> c[j], Val(10)), r.eta) +
        xas_Eqspin_CV(ntuple(j -> c[j + 10], Val(29)), r.eta, S) +
        xas_Uneqspin_CV(ntuple(j -> c[j + 39], Val(10)), r.eta, S, r.chia, r.delta)
end

function xas_amp_cv(i, r, S)
    c = XAS_AMP_FITS[i]
    return xas_Amp_Nospin_CV(ntuple(j -> c[j], Val(8)), r.eta) +
        xas_Amp_Eqspin_CV(ntuple(j -> c[j + 8], Val(28)), r.eta, S) +
        xas_Amp_Uneqspin_CV(ntuple(j -> c[j + 36], Val(6)), r.eta, S, r.chia, r.delta)
end

function xas_inspiral_phase(p, f)
    pn, logs, sigma = p.pn, p.logs, p.sigma
    phi = sum(ntuple(i -> pn[i] * f^((i - 1) / 3), Val(9)))
    phi += logs[1] * f^(5 / 3) * log(f) + logs[2] * f^2 * log(f) + logs[3] * f^(8 / 3) * log(f)
    phi += sum(ntuple(i -> sigma[i] * f^((i + 7) / 3), Val(4)))
    return -(3 * Float64(π)^(-5 / 3) / 128) * phi * f^(-5 / 3)
end

function xas_dinspiral_phase(p, f)
    pn, logs, sigma = p.pn, p.logs, p.sigma
    value = sum(ntuple(i -> pn[i] * ((i - 6) / 3) * f^((i - 9) / 3), Val(9)))
    value += logs[1] / f + logs[2] * f^(-2 / 3) * (1 + log(f) / 3) + logs[3] * (1 + log(f))
    value += sum(ntuple(i -> sigma[i] * ((i + 2) / 3) * f^((i - 1) / 3), Val(4)))
    return -(3 * Float64(π)^(-5 / 3) / 128) * value
end

function xas_rd_phase(p, f)
    c = p.c
    return c[1] * f + 1.5c[2] * f^(2 / 3) - c[3] / f - c[4] / (3f^3) +
        p.cL / p.fdamp * atan((f - p.fRD) / p.fdamp)
end

function xas_drd_phase(p, f)
    c = p.c
    return c[1] + c[2] * f^(-1 / 3) + c[3] / f^2 + c[4] / f^4 + p.cL / (p.fdamp^2 + (f - p.fRD)^2)
end

function xas_int_phase(p, f)
    b = p.b
    return b[1] * f + b[2] * log(f) - b[3] / f - b[4] / (2f^2) - b[5] / (3f^3) +
        2p.cL / p.fdamp * atan((f - p.fRD) / (2p.fdamp))
end

function xas_dint_phase(p, f)
    b = p.b
    return b[1] + b[2] / f + b[3] / f^2 + b[4] / f^3 + b[5] / f^4 + 4p.cL / (4p.fdamp^2 + (f - p.fRD)^2)
end

function xas_phase_coefficients(r, chi1, chi2)
    pn = xas_pn_phase(r.eta, r.delta, chi1, chi2)
    cp = ntuple(i -> 0.0026 + XAS_G4[i] * (1.02r.fMECO - 0.0026), Val(4))
    cv = ntuple(i -> xas_phase_cv(i, r, r.S), Val(4))
    rhs = (cv[1] + cv[3], cv[2] + cv[3], cv[3], cv[4] + cv[3])
    insc = xas_solve(ntuple(i -> (one(cp[i]), cp[i]^(1 / 3), cp[i]^(2 / 3), cp[i]), Val(4)), rhs)
    ins = (; pn = pn[1:9], logs = pn[10:12], sigma = ntuple(i -> -5 / (i + 2) * insc[i], Val(4)))
    fm = 0.6 * (0.5r.fRD + r.fISCO)
    f1 = r.fMECO - 0.03 * (fm - r.fMECO)
    f2 = fm + 0.015 * (fm - r.fMECO)
    rdcp0 = ntuple(i -> fm + XAS_G5[i] * (r.fRD + 1.25r.fdamp - fm), Val(5))
    rdcp = (rdcp0[1], rdcp0[2], rdcp0[3], r.fRD, rdcp0[5])
    v = ntuple(i -> xas_phase_cv(i + 8, r, r.StotR), Val(5))
    rdvalues = (v[1] + v[2] + v[4], v[2] + v[4], v[3] + v[4], v[4], v[5] + v[4])
    rdc = xas_solve(
        ntuple(
            i -> (
                one(rdcp[i]), rdcp[i]^(-1 / 3), rdcp[i]^-2,
                rdcp[i]^-4, -XAS_DPHASE0 / (r.fdamp^2 + (rdcp[i] - r.fRD)^2),
            ), Val(5)
        ), rdvalues
    )
    cL = -XAS_DPHASE0 * rdc[5]
    rd = (; c = rdc[1:4], cL, fRD = r.fRD, fdamp = r.fdamp)
    vi = ntuple(i -> xas_phase_cv(i + 4, r, r.StotR), Val(4))
    intvalues = (
        xas_dinspiral_phase(ins, f1), 0.75 * (vi[1] + v[4]) + 0.25vi[3],
        vi[2] + v[4], vi[4] + vi[2] + v[4], rdvalues[1],
    )
    icp = ntuple(i -> f1 + XAS_G5[i] * (f2 - f1), Val(5))
    irhs = ntuple(i -> intvalues[i] - 4cL / (4r.fdamp^2 + (icp[i] - r.fRD)^2), Val(5))
    bc = xas_solve(ntuple(i -> ntuple(j -> (r.fRD / icp[i])^(j - 1), Val(5)), Val(5)), irhs)
    int = (; b = ntuple(i -> bc[i] * r.fRD^(i - 1), Val(5)), cL, fRD = r.fRD, fdamp = r.fdamp)
    alpha1 = xas_dinspiral_phase(ins, f1) - xas_dint_phase(int, f1)
    alpha0 = xas_inspiral_phase(ins, f1) - xas_int_phase(int, f1) - alpha1 * f1
    beta1 = xas_dint_phase(int, f2) + alpha1 - xas_drd_phase(rd, f2)
    beta0 = xas_int_phase(int, f2) + alpha0 + alpha1 * f2 - xas_rd_phase(rd, f2) - beta1 * f2
    return (; ins, int, rd, f1, f2, alpha0, alpha1, beta0, beta1, eta = r.eta)
end

function phenomxas_phase(p, f)
    value = if f < p.f1
        xas_inspiral_phase(p.ins, f)
    elseif f > p.f2
        xas_rd_phase(p.rd, f) + p.beta0 + p.beta1 * f
    else
        xas_int_phase(p.int, f) + p.alpha0 + p.alpha1 * f
    end
    return value / p.eta
end

function xas_dphase(p, f)
    value = if f < p.f1
        xas_dinspiral_phase(p.ins, f)
    elseif f > p.f2
        xas_drd_phase(p.rd, f) + p.beta1
    else
        xas_dint_phase(p.int, f) + p.alpha1
    end
    return value / p.eta
end

function xas_inspiral_amp(p, f)
    return sum(ntuple(i -> p.pn[i] * f^((i - 1) / 3), Val(7))) +
        sum(ntuple(i -> p.rho[i] * f^((i + 6) / 3), Val(3)))
end

function xas_dinspiral_amp(p, f)
    return sum(ntuple(i -> p.pn[i + 1] * (i / 3) * f^((i - 3) / 3), Val(6))) +
        sum(ntuple(i -> p.rho[i] * ((i + 6) / 3) * f^((i + 3) / 3), Val(3)))
end

xas_rd_amp(p, f) = exp(-(f - p.fRD) * p.gammaR) * p.gammaD13 / ((f - p.fRD)^2 + p.gammaD2)
xas_drd_amp(p, f) = xas_rd_amp(p, f) * (-p.gammaR - 2(f - p.fRD) / ((f - p.fRD)^2 + p.gammaD2))

function xas_amplitude_coefficients(r, chi1, chi2)
    f1 = r.fMECO + 0.25 * (r.fISCO - r.fMECO)
    cv = ntuple(i -> xas_amp_cv(i, r, r.S), Val(3))
    cp = (0.5 * f1, 0.75 * f1, f1)
    rho = xas_solve(ntuple(i -> ntuple(j -> cp[i]^((j + 6) / 3), Val(3)), Val(3)), cv)
    ins = (; pn = xas_pn_amplitude(r.eta, r.delta, chi1, chi2), rho)
    gamma2 = xas_amp_cv(5, r, r.StotR)
    gamma3 = xas_amp_cv(6, r, r.StotR)
    shift = gamma2 <= 1 ? (sqrt(1 - gamma2^2) - 1) / gamma2 : -inv(gamma2)
    f4 = abs(r.fRD + r.fdamp * gamma3 * shift)
    vrd = xas_amp_cv(7, r, r.StotR)
    gamma1 = vrd / (r.fdamp * gamma3) * ((f4 - r.fRD)^2 + (r.fdamp * gamma3)^2) * exp((f4 - r.fRD) * gamma2 / (r.fdamp * gamma3))
    rd = (; fRD = r.fRD, gammaR = gamma2 / (r.fdamp * gamma3), gammaD2 = (gamma3 * r.fdamp)^2, gammaD13 = r.fdamp * gamma1 * gamma3)
    a1, a4 = xas_inspiral_amp(ins, f1), xas_rd_amp(rd, f4)
    v1, v4 = f1^(7 / 6) / a1, f4^(7 / 6) / a4
    v2 = inv(xas_amp_cv(4, r, r.StotR))
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
    m1, m2, chi1, chi2 = order_by_mass(m1in, m2in, chi1in, chi2in)
    r = xas_remnant(m1, m2, chi1, chi2)
    amplitude = xas_amplitude_coefficients(r, chi1, chi2)
    phase_coefficients = xas_phase_coefficients(r, chi1, chi2)
    M_sec = (m1 + m2) * MTSUN_SI
    _, linb, psi4 = xas_phase_at_peak(r.eta, r.StotR, r.chia, r.delta)
    linb -= xas_dphase(phase_coefficients, r.fRD - r.fdamp) + 2Float64(π) * (500 + psi4)
    Mfref = M_sec * fref
    # ChooseFDWaveform's polarization convention includes the additional π.
    phiref = -phenomxas_phase(phase_coefficients, Mfref) - linb * Mfref + 2phase + PI_4 - Float64(π)
    return (; r, amplitude, phase = phase_coefficients, M_sec, linb, phiref)
end

function amplitude_phase_fd(g::IMRPhenomXAS, ::Type{T}, f, p, f_start, f_min, f_max) where {T}
    lo, hi = phenom_frequency_slice(f, f_start, f_min, f_max)
    lo > hi && return (; amplitude = zeros(T, length(f)), phase = zeros(T, length(f)))
    fref = phenom_reference_frequency(g, f_start)
    m = phenomxas_model(p.mass_1, p.mass_2, p.chi_1, p.chi_2, fref, p.phase)
    amp0 = phenom_amplitude_scale(p.mass_1, p.mass_2, p.luminosity_distance)
    norm = amp0 * sqrt(2m.r.eta / 3) * Float64(π)^(-1 / 6)
    cutoff = 0.3 / m.M_sec
    return broadcast_amplitude_phase(phenomxas_signed_amplitude, phenomxas_canonical_phase, T, f, lo, hi, m, norm, cutoff)
end

function phenomxas_signed_amplitude(frequency, model, norm, cutoff)
    frequency <= cutoff || return zero(norm)
    return norm * phenomxas_amplitude(model.amplitude, model.M_sec * frequency)
end

function phenomxas_canonical_phase(frequency, model, amplitude)
    iszero(amplitude) && return zero(amplitude)
    Mf = model.M_sec * frequency
    phi = phenomxas_phase(model.phase, Mf) + model.linb * Mf + model.phiref
    return last(canonical_amplitude_phase(amplitude, phi))
end

# Scalar evaluation at positive frequencies in Hz, using coefficients prepared once per
# waveform. The array assembler handles the configured lower bounds.
function phenomxas_amplitude_phase(frequency, model, norm, cutoff)
    frequency <= cutoff || return (zero(norm), zero(norm))
    Mf = model.M_sec * frequency
    amp = norm * phenomxas_amplitude(model.amplitude, Mf)
    phi = phenomxas_phase(model.phase, Mf) + model.linb * Mf + model.phiref
    return canonical_amplitude_phase(amp, phi)
end

function check_waveform_inputs(g::IMRPhenomXAS, deltaF, p, f_start, f_max)
    (deltaF > 0 && f_start > 0 && f_max > f_start && g.reference_frequency >= 0) || throw(ArgumentError("invalid frequency bounds for IMRPhenomXAS"))
    check_phenom_parameters(p)
    a, b, c, d = order_by_mass(p.mass_1, p.mass_2, p.chi_1, p.chi_2)
    a / b <= 1000 || throw(ArgumentError("IMRPhenomXAS requires mass ratio <= 1000"))
    r = xas_remnant(a, b, c, d)
    (abs(r.a) < 1 && 0 < r.Erad < 1 && r.fMECO <= r.fISCO) || throw(ArgumentError("invalid remnant or MECO/ISCO frequencies"))
    0.3 / ((a + b) * MTSUN_SI) > f_start || throw(ArgumentError("XAS cutoff must exceed the starting frequency"))
    return nothing
end
