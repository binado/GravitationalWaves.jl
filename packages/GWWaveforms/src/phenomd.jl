# IMRPhenomD: aligned-spin inspiral-merger-ringdown model.
#
# Port of LALSimIMRPhenomD.c and LALSimIMRPhenomD_internals.c from LALSimulation (lalsuite
# 7.26, GPL-2.0-or-later); model papers: Husa et al., arXiv:1508.07250 and Khan et al.,
# arXiv:1508.07253. Variable names follow LAL (`fRD`, `fDM`, `delta0`, …) to keep the port
# easy to compare line by line. Frequencies written `Mf` are geometric, `M f`.

# Transition frequencies and cutoff, from LALSimIMRPhenomD.h.
const F_CUT = 0.2
const AMP_FJOIN_INS = 0.014
const PHI_FJOIN_INS = 0.018

# --- LAL's `UsefulPowers` ---------------------------------------------------------------

"""
    UsefulPowers

Powers of a frequency that the model reuses, from a single `pow`: LAL builds them from the
sixth root, and so does this port, to follow its rounding.
"""
struct UsefulPowers{T}
    third::T
    two_thirds::T
    four_thirds::T
    five_thirds::T
    two::T
    seven_thirds::T
    eight_thirds::T
    inv::T
    m_seven_sixths::T
    m_third::T
    m_two_thirds::T
    m_five_thirds::T
end

function useful_powers(number)
    sixth = number^(1.0 / 6.0)
    third = sixth * sixth
    two_thirds = third * third
    four_thirds = number * third
    five_thirds = four_thirds * third
    two = number * number
    seven_thirds = third * two
    eight_thirds = two_thirds * two
    inv = 1.0 / number
    m_sixth = 1.0 / sixth
    m_seven_sixths = inv * m_sixth
    m_third = m_sixth * m_sixth
    m_two_thirds = m_third * m_third
    m_five_thirds = inv * m_two_thirds
    return UsefulPowers(
        third, two_thirds, four_thirds, five_thirds, two, seven_thirds, eight_thirds, inv,
        m_seven_sixths, m_third, m_two_thirds, m_five_thirds,
    )
end

const PI_POWERS = useful_powers(Float64(π))

# --- GSL's natural cubic spline, for the QNM tables -------------------------------------

"""
    NaturalCubicSpline(x, y)

`gsl_interp_cspline`: the natural cubic spline through `(x, y)`, with the coefficients solved
the way GSL does (Cholesky factorisation of the symmetric tridiagonal system). Evaluating it
is plain arithmetic on the knots, so it is differentiable in the evaluation point, and it has
no data-dependent branch or indexing (see [`spline_eval`](@ref)).

`x` holds all knots; `xlo`, `ylo`, `c`, `b` and `d` are the left knot, value and the cubic
coefficients of each of the `n - 1` segments, `y(t) = ylo + δ (b + δ (c + δ d))` with `δ = t - xlo`.
"""
struct NaturalCubicSpline
    x::Vector{Float64}
    xlo::Vector{Float64}
    ylo::Vector{Float64}
    c::Vector{Float64}
    b::Vector{Float64}
    d::Vector{Float64}
end

function NaturalCubicSpline(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    max_index = n - 1
    sys_size = max_index - 1
    c = zeros(n)
    diag = zeros(sys_size)
    offdiag = zeros(sys_size)
    g = zeros(sys_size)
    for i in 1:sys_size
        h_i = x[i + 1] - x[i]
        h_ip1 = x[i + 2] - x[i + 1]
        ydiff_i = y[i + 1] - y[i]
        ydiff_ip1 = y[i + 2] - y[i + 1]
        g_i = h_i != 0.0 ? 1.0 / h_i : 0.0
        g_ip1 = h_ip1 != 0.0 ? 1.0 / h_ip1 : 0.0
        offdiag[i] = h_ip1
        diag[i] = 2.0 * (h_ip1 + h_i)
        g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i)
    end
    # Symmetric tridiagonal solve, A = L D Lᵗ (gsl `solve_tridiag`).
    N = sys_size
    alpha = zeros(N)
    gamma = zeros(N)
    z = zeros(N)
    cc = zeros(N)
    alpha[1] = diag[1]
    gamma[1] = offdiag[1] / alpha[1]
    for i in 2:(N - 1)
        alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1]
        gamma[i] = offdiag[i] / alpha[i]
    end
    alpha[N] = diag[N] - offdiag[N - 1] * gamma[N - 1]
    z[1] = g[1]
    for i in 2:N
        z[i] = g[i] - gamma[i - 1] * z[i - 1]
    end
    for i in 1:N
        cc[i] = z[i] / alpha[i]
    end
    sol = zeros(N)
    sol[N] = cc[N]
    for j in (N - 1):-1:1
        sol[j] = cc[j] - gamma[j] * sol[j + 1]
    end
    c[2:(n - 1)] .= sol
    # The per-segment linear and cubic coefficients, with the expressions GSL evaluates.
    segments = 1:max_index
    dx = x[segments .+ 1] .- x[segments]
    dy = y[segments .+ 1] .- y[segments]
    b = (dy ./ dx) .- dx .* (c[segments .+ 1] .+ 2.0 .* c[segments]) ./ 3.0
    d = (c[segments .+ 1] .- c[segments]) ./ (3.0 .* dx)
    return NaturalCubicSpline(x, x[segments], y[segments], c[segments], b, d)
end

"""
    spline_eval(s::NaturalCubicSpline, x)

Evaluate the spline, extrapolating the first and last cubic outside the knots.
"""
function spline_eval(s::NaturalCubicSpline, x)
    i = clamp(searchsortedlast(s.x, x), 1, length(s.xlo))
    delx = x - s.xlo[i]
    return s.ylo[i] + delx * (s.b[i] + delx * (s.c[i] + delx * s.d[i]))
end

include("phenomd_qnm.jl")
include("phenomd_fits.jl")

const FRING_SPLINE = NaturalCubicSpline(QNM_A, QNM_FRING)
const FDAMP_SPLINE = NaturalCubicSpline(QNM_A, QNM_FDAMP)

# --- final state: arXiv:1508.07250 ------------------------------------------------------

function final_spin_0815_s(eta, s)
    eta2 = eta * eta
    eta3 = eta2 * eta
    s2 = s * s
    s3 = s2 * s
    return eta * (
        3.4641016151377544 - 4.399247300629289 * eta +
            9.397292189321194 * eta2 - 13.180949901606242 * eta3 +
            s * (
            (1.0 / eta - 0.0850917821418767 - 5.837029316602263 * eta) +
                (0.1014665242971878 - 2.0967746996832157 * eta) * s +
                (-1.3546806617824356 + 4.108962025369336 * eta) * s2 +
                (-0.8676969352555539 + 2.064046835273906 * eta) * s3
        )
    )
end

"""
    final_spin_0815(eta, m1, m2, chi1, chi2)

Dimensionless final spin of the remnant (eq. 3.6 of arXiv:1508.07250), for `m1 ≥ m2`.
"""
function final_spin_0815(eta, m1, m2, chi1, chi2)
    M = m1 + m2
    m1s = (m1 / M)^2
    m2s = (m2 / M)^2
    s = (m1s * chi1 + m2s * chi2)
    return final_spin_0815_s(eta, s)
end

function erad_rational_0815_s(eta, s)
    eta2 = eta * eta
    eta3 = eta2 * eta
    return (
        eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta3) *
            (1.0 + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s)
    ) /
        (1.0 + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s)
end

# Fraction of the initial mass radiated (eq. 3.7 and 3.8 of arXiv:1508.07250).
function erad_rational_0815(eta, m1, m2, chi1, chi2)
    M = m1 + m2
    m1s = (m1 / M)^2
    m2s = (m2 / M)^2
    s = (m1s * chi1 + m2s * chi2) / (m1s + m2s)
    return erad_rational_0815_s(eta, s)
end

# Ringdown and damping frequencies, in units of the total initial mass.
fring(erad, finspin) = spline_eval(FRING_SPLINE, finspin) / (1.0 - erad)
fdamp(erad, finspin) = spline_eval(FDAMP_SPLINE, finspin) / (1.0 - erad)

# PN reduced spin parameter, eq. 5.9 of arXiv:1107.1267.
function chi_pn(Seta, eta, chi1, chi2)
    chi_s = (chi1 + chi2)
    chi_a = (chi1 - chi2)
    return 0.5 * (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
end

# --- amplitude --------------------------------------------------------------------------

amp0_func(eta) = sqrt(2.0 / 3.0 * eta) * 0.8263074871107581108331125856317241299 # π^(-1/6)

function amp_ins_ansatz(Mf, pm, pre)
    return 1 + pm.two_thirds * pre.two_thirds +
        pm.four_thirds * pre.four_thirds +
        pm.five_thirds * pre.five_thirds +
        pm.seven_thirds * pre.seven_thirds + pm.eight_thirds * pre.eight_thirds +
        Mf * (pre.one + Mf * pre.two + pm.two * pre.three)
end

function amp_mrd_ansatz(f, a)
    fDMgamma3 = a.fDM * a.gamma3
    fminfRD = f - a.fRD
    return exp(-(fminfRD) * a.gamma2 / (fDMgamma3)) *
        (fDMgamma3 * a.gamma1) / (pow_2_of(fminfRD) + pow_2_of(fDMgamma3))
end

function damp_mrd_ansatz(f, a)
    fDMgamma3 = a.fDM * a.gamma3
    pow2_fDMgamma3 = pow_2_of(fDMgamma3)
    fminfRD = f - a.fRD
    expfactor = exp(((fminfRD) * a.gamma2) / (fDMgamma3))
    pow2pluspow2 = pow_2_of(fminfRD) + pow2_fDMgamma3
    return ((-2 * a.fDM * (fminfRD) * a.gamma3 * a.gamma1) / pow2pluspow2 - (a.gamma2 * a.gamma1)) /
        (expfactor * (pow2pluspow2))
end

# Frequency of the peak of the merger-ringdown amplitude (eq. 20 of arXiv:1508.07253).
function fmax_calc(fRD, fDM, gamma2, gamma3)
    below = !(gamma2 > 1)
    # The square root is only real for `gamma2 ≤ 1`; give it a safe argument otherwise.
    root = sqrt(ifelse(below, 1 - pow_2_of(gamma2), zero(gamma2)))
    return abs(fRD + ifelse(below, (fDM * (-1 + root) * gamma3) / gamma2, (-fDM * gamma3) / gamma2))
end

function amp_int_col_fit_coeff(eta, eta2, chi)
    xi = -1.0 + chi
    return 0.8149838730507785 + 2.5747553517454658 * eta +
        (
        1.1610198035496786 - 2.3627771785551537 * eta + 6.771038707057573 * eta2 +
            (0.7570782938606834 - 2.7256896890432474 * eta + 7.1140380397149965 * eta2) * xi +
            (0.1766934149293479 - 0.7978690983168183 * eta + 2.1162391502005153 * eta2) * xi * xi
    ) * xi
end

function amp_int_ansatz(Mf, a)
    Mf2 = Mf * Mf
    return a.delta0 + Mf * a.delta1 + Mf2 * (a.delta2 + Mf * a.delta3 + a.delta4 * Mf2)
end

# Amplitude coefficients, `ComputeIMRPhenomDAmplitudeCoefficients`; masses are ordered
# `m1 ≥ m2` and `chi1` is the spin on `m1`.
function phenomd_amplitude_coefficients(eta, Seta, chi1, chi2, fRD, fDM)
    etaInv = 1.0 / eta
    chi12 = chi1 * chi1
    chi22 = chi2 * chi2
    eta2 = eta * eta
    eta3 = eta * eta2
    SetaPlus1 = 1.0 + Seta
    chi = chi_pn(Seta, eta, chi1, chi2)
    xi = -1.0 + chi

    gamma1 = gamma1_fun(eta, eta2, xi)
    gamma2 = gamma2_fun(eta, eta2, xi)
    gamma3 = gamma3_fun(eta, eta2, xi)
    fmaxCalc = fmax_calc(fRD, fDM, gamma2, gamma3)
    rho1 = rho1_fun(eta, eta2, xi)
    rho2 = rho2_fun(eta, eta2, xi)
    rho3 = rho3_fun(eta, eta2, xi)

    pre0 = amp_ins_prefactors(eta, chi1, chi2, chi12, chi22, eta2, eta3, Seta, SetaPlus1)
    pre = (; pre0..., seven_thirds = rho1, eight_thirds = rho2, three = rho3)

    # Deltas from three collocation points in [f1, f3].
    f1 = AMP_FJOIN_INS
    f3 = fmaxCalc
    dfx = 0.5 * (f3 - f1)
    f2 = f1 + dfx
    pm1 = useful_powers(f1)
    v1 = amp_ins_ansatz(f1, pm1, pre)
    d1 = damp_ins_ansatz(f1, pm1, eta, chi1, chi2, chi12, chi22, eta2, eta3, Seta, SetaPlus1, rho1, rho2, rho3)
    mrd = (; fRD, fDM, gamma1, gamma2, gamma3)
    v3 = amp_mrd_ansatz(f3, mrd)
    d2 = damp_mrd_ansatz(f3, mrd)
    v2 = amp_int_col_fit_coeff(eta, eta2, chi)

    f12 = f1 * f1
    f13 = f1 * f12
    f14 = f1 * f13
    f15 = f1 * f14
    f22 = f2 * f2
    f23 = f2 * f22
    f24 = f2 * f23
    f32 = f3 * f3
    f33 = f3 * f32
    f34 = f3 * f33
    f35 = f3 * f34
    args = (f1, f2, f3, v1, v2, v3, d1, d2, f12, f13, f14, f15, f22, f23, f24, f32, f33, f34, f35)
    delta0 = delta0_fun(args...)
    delta1 = delta1_fun(args...)
    delta2 = delta2_fun(args...)
    delta3 = delta3_fun(args...)
    delta4 = delta4_fun(args...)

    return (;
        amp0 = amp0_func(eta), pre, fRD, fDM, gamma1, gamma2, gamma3, fmaxCalc,
        delta0, delta1, delta2, delta3, delta4,
    )
end

function _phenomd_amplitude(a, Mf, pm)
    AmpPreFac = a.amp0 * pm.m_seven_sixths
    value = if Mf < AMP_FJOIN_INS
        amp_ins_ansatz(Mf, pm, a.pre)
    elseif Mf >= a.fmaxCalc
        amp_mrd_ansatz(Mf, a)
    else
        amp_int_ansatz(Mf, a)
    end
    return AmpPreFac * value
end

# --- phase ------------------------------------------------------------------------------

function phi_mrd_ansatz_int(f, p)
    sqrootf = sqrt(f)
    fpow1_5 = f * sqrootf
    fpow0_75 = sqrt(fpow1_5)
    return -(p.alpha2 / f) + (4.0 / 3.0) * (p.alpha3 * fpow0_75) + p.alpha1 * f +
        p.alpha4 * atan((f - p.alpha5 * p.fRD) / p.fDM)
end

function dphi_mrd(f, p)
    return (
        p.alpha1 + p.alpha2 / pow_2_of(f) + p.alpha3 / f^0.25 +
            p.alpha4 / (p.fDM * (1 + pow_2_of(f - p.alpha5 * p.fRD) / (pow_2_of(p.fDM))))
    ) * p.etaInv
end

phi_int_ansatz(Mf, p) = p.beta1 * Mf - p.beta3 / (3.0 * pow_3_of(Mf)) + p.beta2 * log(Mf)
dphi_int_ansatz(Mf, p) = (p.beta1 + p.beta3 / pow_4_of(Mf) + p.beta2 / Mf) * p.etaInv
dphi_int_temp(ff, p, C2Int) = C2Int + (p.beta1 + p.beta3 / pow_4_of(ff) + p.beta2 / ff) * p.etaInv

# Prefactors of the inspiral phase that do not depend on frequency (`init_phi_ins_prefactors`).
# `pn.v[k + 1]` is LAL's `pn->v[k]`.
function phi_ins_prefactors(p, pn::TaylorF2Phasing)
    return (;
        initial_phasing = pn.v[6] - PI_4,
        two_thirds = pn.v[8] * PI_POWERS.two_thirds,
        third = pn.v[7] * PI_POWERS.third,
        third_with_logv = pn.vlogv6 * PI_POWERS.third,
        logv = pn.vlogv5,
        minus_third = pn.v[5] * PI_POWERS.m_third,
        minus_two_thirds = pn.v[4] * PI_POWERS.m_two_thirds,
        minus_one = pn.v[3] * PI_POWERS.inv,
        minus_four_thirds = pn.v[2] / PI_POWERS.four_thirds,
        minus_five_thirds = pn.v[1] * PI_POWERS.m_five_thirds,
        one = p.sigma1,
        four_thirds = p.sigma2 * 0.75,
        five_thirds = p.sigma3 * 0.6,
        two = p.sigma4 * 0.5,
    )
end

function phi_ins_ansatz_int(Mf, pm, pre, p)
    v = pm.third * PI_POWERS.third
    logv = log(v)
    phasing = pre.initial_phasing
    phasing += pre.two_thirds * pm.two_thirds
    phasing += pre.third * pm.third
    phasing += pre.third_with_logv * logv * pm.third
    phasing += pre.logv * logv
    phasing += pre.minus_third * pm.m_third
    phasing += pre.minus_two_thirds * pm.m_two_thirds
    phasing += pre.minus_one * pm.inv
    phasing += pre.minus_four_thirds / pm.four_thirds
    phasing += pre.minus_five_thirds * pm.m_five_thirds
    # Higher-order terms calibrated for PhenomD.
    phasing += (
        pre.one * Mf + pre.four_thirds * pm.four_thirds + pre.five_thirds * pm.five_thirds +
            pre.two * pm.two
    ) * p.etaInv
    return phasing
end

function dphi_ins_ansatz_int(Mf, p, pn::TaylorF2Phasing)
    Pi = Float64(π)
    v = cbrt(Pi * Mf)
    logv = log(v)
    v2 = v * v
    v3 = v * v2
    v4 = v * v3
    v5 = v * v4
    v6 = v * v5
    v7 = v * v6
    v8 = v * v7
    Dphasing = zero(v)
    Dphasing += 2.0 * pn.v[8] * v7
    Dphasing += (pn.v[7] + pn.vlogv6 * (1.0 + logv)) * v6
    Dphasing += pn.vlogv5 * v5
    Dphasing += -1.0 * pn.v[5] * v4
    Dphasing += -2.0 * pn.v[4] * v3
    Dphasing += -3.0 * pn.v[3] * v2
    Dphasing += -4.0 * pn.v[2] * v
    Dphasing += -5.0 * pn.v[1]
    Dphasing /= v8 * 3.0
    Dphasing *= Pi
    Dphasing += (
        p.sigma1 + p.sigma2 * v * PI_POWERS.m_third + p.sigma3 * v2 * PI_POWERS.m_two_thirds +
            (p.sigma4 * PI_POWERS.inv) * v3
    ) * p.etaInv
    return Dphasing
end

# Phase coefficients before the C¹ continuity constants (`ComputeIMRPhenomDPhaseCoefficients`).
function phenomd_phase_fits(eta, Seta, chi1, chi2, fRD, fDM)
    eta2 = eta * eta
    chi = chi_pn(Seta, eta, chi1, chi2)
    xi = -1.0 + chi
    return (;
        etaInv = 1.0 / eta,
        sigma1 = sigma1Fit(eta, eta2, xi), sigma2 = sigma2Fit(eta, eta2, xi),
        sigma3 = sigma3Fit(eta, eta2, xi), sigma4 = sigma4Fit(eta, eta2, xi),
        beta1 = beta1Fit(eta, eta2, xi), beta2 = beta2Fit(eta, eta2, xi), beta3 = beta3Fit(eta, eta2, xi),
        alpha1 = alpha1Fit(eta, eta2, xi), alpha2 = alpha2Fit(eta, eta2, xi),
        alpha3 = alpha3Fit(eta, eta2, xi), alpha4 = alpha4Fit(eta, eta2, xi),
        alpha5 = alpha5Fit(eta, eta2, xi),
        fRD, fDM,
    )
end

# Join the inspiral, intermediate and merger-ringdown phases with continuous value and
# first derivative (`ComputeIMRPhenDPhaseConnectionCoefficients`).
function phase_connection(p, pn, pre)
    etaInv = p.etaInv
    fMRDJoin = 0.5 * p.fRD

    DPhiIns = dphi_ins_ansatz_int(PHI_FJOIN_INS, p, pn)
    DPhiInt = dphi_int_ansatz(PHI_FJOIN_INS, p)
    C2Int = DPhiIns - DPhiInt
    pm = useful_powers(PHI_FJOIN_INS)
    C1Int = phi_ins_ansatz_int(PHI_FJOIN_INS, pm, pre, p) -
        etaInv * phi_int_ansatz(PHI_FJOIN_INS, p) - C2Int * PHI_FJOIN_INS

    PhiIntTempVal = etaInv * phi_int_ansatz(fMRDJoin, p) + C1Int + C2Int * fMRDJoin
    DPhiIntTempVal = dphi_int_temp(fMRDJoin, p, C2Int)
    DPhiMRDVal = dphi_mrd(fMRDJoin, p)
    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - etaInv * phi_mrd_ansatz_int(fMRDJoin, p) - C2MRD * fMRDJoin
    return (; C1Int, C2Int, C1MRD, C2MRD, fMRDJoin)
end

function _phenomd_phase(m, Mf, pm)
    p = m.phi
    if Mf < PHI_FJOIN_INS
        return phi_ins_ansatz_int(Mf, pm, m.phi_pre, p)
    elseif Mf >= m.conn.fMRDJoin
        return p.etaInv * phi_mrd_ansatz_int(Mf, p) + m.conn.C1MRD + m.conn.C2MRD * Mf
    else
        return p.etaInv * phi_int_ansatz(Mf, p) + m.conn.C1Int + m.conn.C2Int * Mf
    end
end

# 3PN spin-spin term that LAL's TaylorF2 has but PhenomD was tuned without.
function subtract_3pn_ss(m1, m2, M, eta, chi1, chi2)
    m1M = m1 / M
    m2M = m2 / M
    pn_ss3 = (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1 * chi2
    pn_ss3 += ((4703.5 / 8.4 + 2935.0 / 6.0 * m1M - 120.0 * m1M * m1M) + (-4108.25 / 6.72 - 108.5 / 1.2 * m1M + 125.5 / 3.6 * m1M * m1M)) * m1M * m1M * chi1 * chi1
    pn_ss3 += ((4703.5 / 8.4 + 2935.0 / 6.0 * m2M - 120.0 * m2M * m2M) + (-4108.25 / 6.72 - 108.5 / 1.2 * m2M + 125.5 / 3.6 * m2M * m2M)) * m2M * m2M * chi2 * chi2
    return pn_ss3
end

"""
    PhenomDModel

Everything the IMRPhenomD strain needs that does not depend on the frequency, built once per
waveform by [`phenomd_model`](@ref).
"""
struct PhenomDModel{A, P, N, PF, C, T}
    M_sec::T
    amp::A
    phi::P
    pn::N
    phi_pre::PF
    conn::C
    t0::T
    phi_precalc::T
    MfRef::T
end

# Symmetric mass ratio with `PhenomInternal_nudge`: eta is snapped to 1/4 if it exceeds it only
# through rounding. Larger values are rejected by `check_waveform_inputs`.
function nudged_eta(m1, m2)
    M = m1 + m2
    eta = m1 * m2 / (M * M)
    return ifelse(eta > 0.25, oftype(eta, 0.25), eta)
end

"""
    phenomd_model(m1, m2, chi1, chi2, f_ref, phi0) -> PhenomDModel

Set up IMRPhenomD for component masses `m1`, `m2` (solar masses, any order) with aligned
spins `chi1`, `chi2`, fixing the orbital phase `phi0` at the reference frequency `f_ref`
(Hz). Mirrors `IMRPhenomDGenerateFD`, including the swap to `m1 ≥ m2`. The inputs are trusted and
[`check_waveform_inputs`](@ref) checks eta and the final spin before numerical evaluation.
"""
function phenomd_model(m1_in, m2_in, chi1_in, chi2_in, fRef, phi0)
    m1, m2, chi1, chi2 = order_by_mass(m1_in, m2_in, chi1_in, chi2_in)
    M = m1 + m2
    eta = nudged_eta(m1, m2)
    M_sec = M * MTSUN_SI

    # Carry the mass difference directly: sqrt(1 - 4eta) has singular AD at equal mass.
    Seta = (m1 - m2) / M
    finspin = final_spin_0815(eta, m1, m2, chi1, chi2)
    erad = erad_rational_0815(eta, m1, m2, chi1, chi2)
    fRD = fring(erad, finspin)
    fDM = fdamp(erad, finspin)

    amp = phenomd_amplitude_coefficients(eta, Seta, chi1, chi2, fRD, fDM)
    phi = phenomd_phase_fits(eta, Seta, chi1, chi2, fRD, fDM)

    # TaylorF2 phasing to 3.5PN spin order, minus the 3PN spin-spin term (see `subtract_3pn_ss`).
    pn0 = taylorf2_pn_phasing(m1, m2, chi1, chi2; spin_order = 7)
    testGRcor = 1.0
    v6 = pn0.v[7] - (subtract_3pn_ss(m1, m2, M, eta, chi1, chi2) * pn0.v[1]) * testGRcor
    pn = TaylorF2Phasing(
        (pn0.v[1], pn0.v[2], pn0.v[3], pn0.v[4], pn0.v[5], pn0.v[6], v6, pn0.v[8]),
        pn0.vlogv5, pn0.vlogv6,
    )

    phi_pre = phi_ins_prefactors(phi, pn)
    conn = phase_connection(phi, pn, phi_pre)
    phi = (; phi..., conn...)

    # Time shift so that the peak amplitude sits near t = 0.
    t0 = dphi_mrd(amp.fmaxCalc, phi)

    MfRef = M_sec * fRef
    pmRef = useful_powers(MfRef)
    m0 = (; phi, phi_pre, conn)
    phifRef = _phenomd_phase(m0, MfRef, pmRef)
    # Factor of 2 because phi0 is the orbital phase.
    phi_precalc = 2.0 * phi0 + phifRef
    return PhenomDModel(M_sec, amp, phi, pn, phi_pre, conn, t0, phi_precalc, MfRef)
end

"""
    phenomd_amplitude(model::PhenomDModel, Mf)

IMRPhenomD amplitude at the geometric frequency `Mf` (the factor multiplying `amp0` in
`h̃ = amp0 · A(Mf) e^{-iΨ(Mf)}`).
"""
phenomd_amplitude(model::PhenomDModel, Mf) = _phenomd_amplitude(model.amp, Mf, useful_powers(Mf))

"""
    phenomd_phase(model::PhenomDModel, Mf)

IMRPhenomD phase Ψ(Mf) before the alignment at the reference frequency.
"""
phenomd_phase(model::PhenomDModel, Mf) = _phenomd_phase(model, Mf, useful_powers(Mf))

# --- driver -----------------------------------------------------------------------------

"""
    phenomd_fd(T, nout, deltaF, p, g, f_start, f_max)

Frequency-domain IMRPhenomD strain for the optimally oriented source on `f_i = i deltaF`,
`i = 0 … nout - 1`, as `XLALSimIMRPhenomDGenerateFD`: nonzero from `f_start` up to
`min(f_max, 0.2 / M)` (and `f_max == 0` means up to `0.2 / M`), and the sample at the upper
end itself is zero, as in LAL. The inclination factors are applied by `hphc_fd`.

The inputs are trusted (see [`validate_parameters`](@ref)).
This is the reference loop; [`hphc_fd`](@ref) uses [`phenomd_amplitude_phase_slice`](@ref) to apply the public frequency bounds.
"""
function phenomd_fd(::Type{T}, nout::Integer, deltaF, p, g, f_start, f_max) where {T}
    m1 = p.mass_1
    m2 = p.mass_2
    fRef = phenom_reference_frequency(g, f_start)
    fCut = F_CUT / ((m1 + m2) * MTSUN_SI)
    f_max_prime = f_max == 0 ? fCut : f_max
    f_max_prime = f_max_prime > fCut ? fCut : f_max_prime

    model = phenomd_model(m1, m2, p.chi_1, p.chi_2, fRef, p.phase)
    amp0 = phenom_amplitude_scale(m1, m2, p.luminosity_distance)

    i_start = trunc(Int, f_start / deltaF)
    i_stop = trunc(Int, f_max_prime / deltaF)
    h = zeros(Complex{T}, nout)
    for i in i_start:(min(i_stop, nout) - 1)
        Mf = model.M_sec * (i * deltaF)
        pm = useful_powers(Mf)
        amp = _phenomd_amplitude(model.amp, Mf, pm)
        phi = _phenomd_phase(model, Mf, pm)
        phi -= model.t0 * (Mf - model.MfRef) + model.phi_precalc
        a = amp0 * amp
        h[i + 1] = complex(a * cos(phi), -(a * sin(phi)))
    end
    return h
end

"""
    phenomd_amplitude_phase_slice(T, f, lo, hi, p, g, f_start, f_max)

Nonnegative amplitude and unwrapped phase of [`phenomd_fd`](@ref), restricted to the one-based indices `lo:hi` from
[`phenom_frequency_slice`](@ref). Frequencies are read directly from `f`. The upper edge
is rounded down to a multiple of `step(f)` and excluded, preserving LAL's cutoff convention.
"""
function phenomd_amplitude_phase_slice(
        ::Type{T}, f::AbstractRange, lo::Integer, hi::Integer, p, g, f_start, f_max
    ) where {T}
    lo > hi && return (; amplitude = zeros(T, length(f)), phase = zeros(T, length(f)))
    m1 = p.mass_1
    m2 = p.mass_2
    fRef = phenom_reference_frequency(g, f_start)
    M = m1 + m2
    fCut = F_CUT / (M * MTSUN_SI)
    f_max_eff = ifelse(fCut < f_max, fCut, f_max)

    model = phenomd_model(m1, m2, p.chi_1, p.chi_2, fRef, p.phase)
    amp0 = phenom_amplitude_scale(m1, m2, p.luminosity_distance)

    deltaF = Float64(step(f))
    stop = floor(f_max_eff / deltaF) * deltaF
    return broadcast_amplitude_phase(phenomd_signed_amplitude, phenomd_canonical_phase, T, f, lo, hi, model, amp0, stop)
end

# Separate scalar outputs allow broadcasting directly into the public arrays. The phase
# receives the signed amplitude so canonicalization does not evaluate the amplitude again.
function phenomd_signed_amplitude(frequency, model, amp0, cutoff)
    frequency < cutoff || return zero(amp0)
    Mf = model.M_sec * frequency
    return amp0 * _phenomd_amplitude(model.amp, Mf, useful_powers(Mf))
end

function phenomd_canonical_phase(frequency, model, amplitude)
    iszero(amplitude) && return zero(amplitude)
    Mf = model.M_sec * frequency
    phi = _phenomd_phase(model, Mf, useful_powers(Mf)) - (model.t0 * (Mf - model.MfRef) + model.phi_precalc)
    return last(canonical_amplitude_phase(amplitude, -phi))
end

# Scalar evaluation at positive frequencies in Hz. Coefficients and the LAL-rounded upper
# cutoff are prepared once per waveform; the array assembler handles the lower bounds.
function phenomd_amplitude_phase(frequency, model, amp0, cutoff)
    frequency < cutoff || return (zero(amp0), zero(amp0))
    Mf = model.M_sec * frequency
    pm = useful_powers(Mf)
    amp = amp0 * _phenomd_amplitude(model.amp, Mf, pm)
    phi = _phenomd_phase(model, Mf, pm) - (model.t0 * (Mf - model.MfRef) + model.phi_precalc)
    return canonical_amplitude_phase(amp, -phi)
end

# Inputs `phenomd_fd` does not check; run by `validate_parameters`.
function check_waveform_inputs(g::IMRPhenomD, deltaF, p, f_start, f_max)
    m1 = p.mass_1
    m2 = p.mass_2
    fRef_in = g.reference_frequency
    fRef_in >= 0 || throw(ArgumentError("reference frequency must not be negative"))
    deltaF > 0 || throw(ArgumentError("frequency step must be positive"))
    f_start > 0 || throw(ArgumentError("minimum frequency must be positive"))
    f_max > 0 || throw(ArgumentError(MAX_FREQUENCY_ISCO_MESSAGE))
    check_phenom_parameters(p)

    fCut = F_CUT / ((m1 + m2) * MTSUN_SI)
    fCut > f_start || throw(ArgumentError("fCut = $fCut Hz is not above f_min = $f_start Hz"))
    f_max_prime = f_max > fCut ? fCut : f_max
    f_max_prime > f_start || throw(ArgumentError("f_max <= f_min"))

    a, b, chi1, chi2 = order_by_mass(m1, m2, p.chi_1, p.chi_2)
    eta = a * b / ((a + b) * (a + b))
    eta - 0.25 <= 5.0e-7 || throw(ArgumentError("unphysical eta; must be between 0 and 0.25"))
    final_spin_0815(nudged_eta(a, b), a, b, chi1, chi2) > 1.0 &&
        throw(ArgumentError("final spin > 1.0 is not supported"))
    return nothing
end
