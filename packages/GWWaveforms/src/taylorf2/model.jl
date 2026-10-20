# TaylorF2: stationary-phase post-Newtonian inspiral for aligned spins.
#
# Port of `XLALSimInspiralPNPhasing_F2` (LALSimInspiralPNCoefficients.c) and
# `XLALSimInspiralTaylorF2Core` / `XLALSimInspiralTaylorF2` (LALSimInspiralTaylorF2.c) from
# LALSimulation (lalsuite 7.26), GPL-2.0-or-later. Tidal terms are dropped (the tidal
# deformabilities are zero) and so are the testing-GR deviations of the phasing coefficients.

"""
    TaylorF2Phasing(v, vlogv5, vlogv6)

Coefficients of the SPA phase `Σ_k v[k+1] v^k + vlogv_k v^k log v`, `k = 0…7`, in powers of
`v = (π M f)^{1/3}`, including the overall `3 / (128 η)`. The phase is `phasing / v^5`.
"""
struct TaylorF2Phasing{T}
    v::NTuple{8, T}
    vlogv5::T
    vlogv6::T
end

# --- coefficients, as in LALSimInspiralPNCoefficients.c ---------------------------------

taylorf2_phasing_2pn(eta) = 5 * (74.3 / 8.4 + 11 * eta) / 9
taylorf2_phasing_3pn(eta) = -16 * Float64(π)
taylorf2_phasing_4pn(eta) = 5 * (3058.673 / 7.056 + 5429 / 7 * eta + 617 * eta * eta) / 72
taylorf2_phasing_5pn(eta) = 5 / 9 * (772.9 / 8.4 - 13 * eta) * Float64(π)
taylorf2_phasing_5pn_log(eta) = 5 / 3 * (772.9 / 8.4 - 13 * eta) * Float64(π)
taylorf2_phasing_6pn_log(eta) = -684.8 / 2.1
function taylorf2_phasing_6pn(eta)
    return 11583.231236531 / 4.69421568 - 640 / 3 * Float64(π) * Float64(π) -
        684.8 / 2.1 * GAMMA + eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * Float64(π) * Float64(π)) +
        eta * eta * 76.055 / 1.728 - eta * eta * eta * 127.825 / 1.296 +
        taylorf2_phasing_6pn_log(eta) * log(4.0)
end
function taylorf2_phasing_7pn(eta)
    return Float64(π) * (770.96675 / 2.54016 + 378.515 / 1.512 * eta - 740.45 / 7.56 * eta * eta)
end

# Spin-orbit terms, with `m` the mass of the body over the total mass.
taylorf2_phasing_3pn_so(m) = m * (25 + 38 / 3 * m)
function taylorf2_phasing_5pn_so(m)
    return -m * (1391.5 / 8.4 - m * (1 - m) * 10 / 3 + m * (1276 / 8.1 + m * (1 - m) * 170 / 9))
end
taylorf2_phasing_6pn_so(m) = Float64(π) * m * (1490 / 3 + m * 260)
function taylorf2_phasing_7pn_so(m)
    eta = m * (1 - m)
    return m * (
        -17097.8035 / 4.8384 + eta * 28764.25 / 6.72 + eta * eta * 47.35 / 1.44 +
            m * (-7189.233785 / 1.524096 + eta * 458.555 / 3.024 - eta * eta * 534.5 / 7.2)
    )
end

# Spin-spin, quadrupole-monopole and self-spin terms.
taylorf2_phasing_4pn_s1s2(eta) = 247 / 4.8 * eta
taylorf2_phasing_4pn_s1s2o(eta) = -721 / 4.8 * eta
taylorf2_phasing_4pn_qm2so(m) = -720 / 9.6 * m * m
taylorf2_phasing_4pn_self2so(m) = 1 / 9.6 * m * m
taylorf2_phasing_4pn_qm2s(m) = 240 / 9.6 * m * m
taylorf2_phasing_4pn_self2s(m) = -7 / 9.6 * m * m
taylorf2_phasing_6pn_s1s2o(eta) = (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta
taylorf2_phasing_6pn_self2s(m) = (-4108.25 / 6.72 - 108.5 / 1.2 * m + 125.5 / 3.6 * m * m) * m * m
taylorf2_phasing_6pn_qm2s(m) = (4703.5 / 8.4 + 2935 / 6 * m - 120 * m * m) * m * m

"""
    taylorf2_pn_phasing(m1, m2, chi1, chi2; spin_order = -1) -> TaylorF2Phasing

Port of `XLALSimInspiralPNPhasing_F2` for black holes (no tides, quadrupole-monopole
parameter 1) with aligned dimensionless spins `chi1`, `chi2` and masses in solar masses.
`spin_order` is LAL's `LAL_SIM_INSPIRAL_SPIN_ORDER` in half post-Newtonian orders:
`-1` or `7` for everything up to 3.5PN, down to `0` for no spin terms.
"""
function taylorf2_pn_phasing(m1, m2, chi1, chi2; spin_order::Integer = -1)
    chi1sq = chi1 * chi1
    chi2sq = chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    mtot = m1 + m2
    eta = m1 * m2 / mtot / mtot
    m1_m = m1 / mtot
    m2_m = m2 / mtot
    pfa_n = 3 / (128 * eta)

    v0 = one(eta)
    v1 = zero(eta)
    v2 = taylorf2_phasing_2pn(eta)
    v3 = taylorf2_phasing_3pn(eta)
    v4 = taylorf2_phasing_4pn(eta)
    v5 = taylorf2_phasing_5pn(eta)
    l5 = taylorf2_phasing_5pn_log(eta)
    v6 = taylorf2_phasing_6pn(eta)
    l6 = taylorf2_phasing_6pn_log(eta)
    v7 = taylorf2_phasing_7pn(eta)

    # The spin orders cascade: each includes all the lower ones.
    so = spin_order == -1 ? 7 : spin_order
    qm1 = 1.0 # 1 + d_quad_mon, zero for black holes
    qm2 = 1.0
    if so >= 7
        v7 += taylorf2_phasing_7pn_so(m1_m) * chi1 + taylorf2_phasing_7pn_so(m2_m) * chi2
    end
    if so >= 6
        v6 += taylorf2_phasing_6pn_so(m1_m) * chi1 + taylorf2_phasing_6pn_so(m2_m) * chi2 +
            taylorf2_phasing_6pn_s1s2o(eta) * chi1 * chi2 +
            (taylorf2_phasing_6pn_qm2s(m1_m) * qm1 + taylorf2_phasing_6pn_self2s(m1_m)) * chi1sq +
            (taylorf2_phasing_6pn_qm2s(m2_m) * qm2 + taylorf2_phasing_6pn_self2s(m2_m)) * chi2sq
    end
    if so >= 5
        so5 = taylorf2_phasing_5pn_so(m1_m) * chi1 + taylorf2_phasing_5pn_so(m2_m) * chi2
        v5 += so5
        l5 += 3 * so5
    end
    if so >= 4
        v4 += taylorf2_phasing_4pn_s1s2(eta) * chi1dotchi2 + taylorf2_phasing_4pn_s1s2o(eta) * chi1 * chi2 +
            (taylorf2_phasing_4pn_qm2so(m1_m) * qm1 + taylorf2_phasing_4pn_self2so(m1_m)) * chi1 * chi1 +
            (taylorf2_phasing_4pn_qm2so(m2_m) * qm2 + taylorf2_phasing_4pn_self2so(m2_m)) * chi2 * chi2 +
            (taylorf2_phasing_4pn_qm2s(m1_m) * qm1 + taylorf2_phasing_4pn_self2s(m1_m)) * chi1sq +
            (taylorf2_phasing_4pn_qm2s(m2_m) * qm2 + taylorf2_phasing_4pn_self2s(m2_m)) * chi2sq
    end
    if so >= 3
        v3 += taylorf2_phasing_3pn_so(m1_m) * chi1 + taylorf2_phasing_3pn_so(m2_m) * chi2
    end

    return TaylorF2Phasing(
        (v0 * pfa_n, v1 * pfa_n, v2 * pfa_n, v3 * pfa_n, v4 * pfa_n, v5 * pfa_n, v6 * pfa_n, v7 * pfa_n),
        l5 * pfa_n, l6 * pfa_n,
    )
end

"""
    truncate_phasing(pn::TaylorF2Phasing, phase_order) -> TaylorF2Phasing

Drop the terms above `phase_order` (in whole post-Newtonian orders, `-1` or `7` keep all), as
`XLALSimInspiralTaylorF2Core` does.
"""
function truncate_phasing(pn::TaylorF2Phasing{T}, phase_order::Integer) where {T}
    po = phase_order == -1 ? 7 : phase_order
    z = zero(T)
    v = ntuple(k -> k - 1 <= po ? pn.v[k] : z, Val(8))
    return TaylorF2Phasing(v, po >= 5 ? pn.vlogv5 : z, po >= 6 ? pn.vlogv6 : z)
end

"""
    taylorf2_phasing(pn::TaylorF2Phasing, v)

The SPA phase `Ψ(f)` without time and phase shifts, at `v = (π M f)^{1/3}`.
"""
function taylorf2_phasing(pn::TaylorF2Phasing, v)
    logv = log(v)
    v2 = v * v
    v3 = v * v2
    v4 = v * v3
    v5 = v * v4
    v6 = v * v5
    v7 = v * v6
    phasing = zero(v)
    phasing += pn.v[8] * v7
    phasing += (pn.v[7] + pn.vlogv6 * logv) * v6
    phasing += (pn.v[6] + pn.vlogv5 * logv) * v5
    phasing += pn.v[5] * v4
    phasing += pn.v[4] * v3
    phasing += pn.v[3] * v2
    phasing += pn.v[2] * v
    phasing += pn.v[1]
    return phasing / v5
end

# --- amplitude, SPA corrections to the Newtonian amplitude ------------------------------

taylorf2_energy_0pn(eta) = -eta / 2.0
taylorf2_energy_2pn(eta) = -(0.75 + eta / 12.0)
taylorf2_energy_4pn(eta) = -(27.0 / 8.0 - 19.0 / 8.0 * eta + 1.0 / 24.0 * eta * eta)
function taylorf2_energy_6pn(eta)
    return -(
        67.5 / 6.4 - (344.45 / 5.76 - 20.5 / 9.6 * Float64(π) * Float64(π)) * eta +
            15.5 / 9.6 * eta * eta + 3.5 / 518.4 * eta * eta * eta
    )
end

taylorf2_flux_0pn(eta) = 32.0 * eta * eta / 5.0
taylorf2_flux_2pn(eta) = -(12.47 / 3.36 + 3.5 / 1.2 * eta)
taylorf2_flux_3pn(eta) = 4.0 * Float64(π)
taylorf2_flux_4pn(eta) = -(44.711 / 9.072 - 92.71 / 5.04 * eta - 6.5 / 1.8 * eta * eta)
taylorf2_flux_5pn(eta) = -(81.91 / 6.72 + 58.3 / 2.4 * eta) * Float64(π)
function taylorf2_flux_6pn(eta)
    return (
        664.3739519 / 6.98544 + 16.0 / 3.0 * Float64(π) * Float64(π) - 17.12 / 1.05 * GAMMA -
            17.12 / 1.05 * log(4.0) + (4.1 / 4.8 * Float64(π) * Float64(π) - 134.543 / 7.776) * eta -
            94.403 / 3.024 * eta * eta - 7.75 / 3.24 * eta * eta * eta
    )
end
taylorf2_flux_6pn_log(eta) = -17.12 / 1.05
taylorf2_flux_7pn(eta) = -(162.85 / 5.04 - 214.745 / 1.728 * eta - 193.385 / 3.024 * eta * eta) * Float64(π)

"""
    TaylorF2Amplitude(eta, amplitude_order)

Energy and flux coefficients of the SPA amplitude corrections, `amplitude_order` in
`{-1, 0, 2, …, 7}` (`-1` is the same as `0`: no corrections by default).
"""
struct TaylorF2Amplitude{T}
    order::Int
    f_ta_n::T
    FTa2::T
    FTa3::T
    FTa4::T
    FTa5::T
    FTl6::T
    FTa6::T
    FTa7::T
    d_e_ta_n::T
    d_e_ta1::T
    d_e_ta2::T
    d_e_ta3::T
end

function TaylorF2Amplitude(eta, amplitude_order::Integer)
    # Some coefficients are pure numbers; promote so all share the type of `eta`.
    c = promote(
        taylorf2_flux_0pn(eta), taylorf2_flux_2pn(eta), taylorf2_flux_3pn(eta),
        taylorf2_flux_4pn(eta), taylorf2_flux_5pn(eta), taylorf2_flux_6pn_log(eta),
        taylorf2_flux_6pn(eta), taylorf2_flux_7pn(eta),
        2.0 * taylorf2_energy_0pn(eta), 2.0 * taylorf2_energy_2pn(eta),
        3.0 * taylorf2_energy_4pn(eta), 4.0 * taylorf2_energy_6pn(eta),
    )
    return TaylorF2Amplitude(Int(amplitude_order), c...)
end

# sqrt(-d_e/F) with the corrections up to `order`; multiplies the Newtonian amplitude.
function taylorf2_amplitude_factor(a::TaylorF2Amplitude, v)
    logv = log(v)
    v2 = v * v
    v3 = v * v2
    v4 = v * v3
    v5 = v * v4
    v6 = v * v5
    v7 = v * v6
    v8 = v * v7
    v9 = v * v8
    v10 = v * v9
    flux = zero(v)
    d_energy = zero(v)
    o = a.order
    if o >= 7
        flux += a.FTa7 * v7
    end
    if o >= 6
        flux += (a.FTa6 + a.FTl6 * logv) * v6
        d_energy += a.d_e_ta3 * v6
    end
    if o >= 5
        flux += a.FTa5 * v5
    end
    if o >= 4
        flux += a.FTa4 * v4
        d_energy += a.d_e_ta2 * v4
    end
    if o >= 3
        flux += a.FTa3 * v3
    end
    if o >= 2
        flux += a.FTa2 * v2
        d_energy += a.d_e_ta1 * v2
    end
    flux += 1.0
    d_energy += 1.0
    flux *= a.f_ta_n * v10
    d_energy *= a.d_e_ta_n * v
    return sqrt(-d_energy / flux)
end

# --- driver ------------------------------------------------------------------------------

# LAL stores the coalescence time shift -1/delta_f in a `LIGOTimeGPS`, whole seconds plus
# nanoseconds rounded to the nearest integer, and builds the phase shift from that.
function gps_time_shift(delta_f)
    t = -1 / delta_f
    ns = round(Int64, rem(t, 1.0) * 1.0e9)
    total = Int64(trunc(t)) * 1_000_000_000 + ns
    seconds = fld(total, 1_000_000_000)
    nanoseconds = mod(total, 1_000_000_000)
    return 2 * Float64(π) * (seconds + 1.0e-9 * nanoseconds)
end

# Inputs `amplitude_phase_fd` does not check; run by `validate_parameters`.
function check_waveform_inputs(g::TaylorF2, delta_f, params::NamedTuple, f_start, f_max)
    m1 = params.mass_1
    m2 = params.mass_2
    (m1 > 0 && m2 > 0) || throw(ArgumentError("masses must be positive"))
    f_start > 0 || throw(ArgumentError("minimum frequency must be positive"))
    g.reference_frequency >= 0 || throw(ArgumentError("reference frequency must not be negative"))
    params.luminosity_distance > 0 || throw(ArgumentError("luminosity distance must be positive"))
    delta_f > 0 || throw(ArgumentError("frequency step must be positive"))
    f_max > 0 || throw(ArgumentError(MAX_FREQUENCY_ISCO_MESSAGE))
    f_max > f_start || throw(
        ArgumentError("maximum frequency $f_max must exceed the start frequency $f_start")
    )
    return nothing
end

# TaylorF2 (`XLALSimInspiralTaylorF2`) on the supplied frequencies: nonnegative amplitude and
# unwrapped phase, zero outside `[f_start, f_max]`. Grids may start above zero; the loop never
# evaluates the model at `f = 0`. The inclination factors are applied by `hphc_fd`. Inputs are
# trusted; run `validate_parameters` on the host first.
function amplitude_phase_fd(f::AbstractRange{Float64}, params::NamedTuple, g::TaylorF2)
    T, p = promote_parameters(aligned_parameters(params))
    f_min, f_max = frequency_bounds(g, f)
    f_start = start_frequency(g, f_min)
    lo, hi = frequency_slice(f, f_start, f_min, f_max)
    delta_f = step(f)
    amplitude = zeros(T, length(f))
    phase = zeros(T, length(f))
    lo > hi && return (; amplitude, phase)
    m1 = p.mass_1
    m2 = p.mass_2
    r = p.luminosity_distance * 1.0e6 * PC_SI
    f_ref = g.reference_frequency

    m = m1 + m2
    m_sec = m * MTSUN_SI
    eta = m1 * m2 / (m * m)
    pi_m = Float64(π) * m_sec
    shft = gps_time_shift(delta_f)

    pn = truncate_phasing(
        taylorf2_pn_phasing(m1, m2, p.chi_1, p.chi_2; spin_order = g.pn_spin_order),
        g.pn_phase_order,
    )
    correction = TaylorF2Amplitude(eta, g.pn_amplitude_order)
    amp0 = -4.0 * m1 * m2 / r * MRSUN_SI * MTSUN_SI * sqrt(Float64(π) / 12.0)
    ref_phasing = f_ref != 0 ? taylorf2_phasing(pn, cbrt(pi_m * f_ref)) : zero(amp0)

    for i in lo:hi
        frequency = f[i]
        v = cbrt(pi_m * frequency)
        phasing = taylorf2_phasing(pn, v)
        amp = amp0 * taylorf2_amplitude_factor(correction, v) * v
        x = (phasing + ((shft * frequency - 2.0 * p.phase) - ref_phasing)) - PI_4
        amplitude[i], phase[i] = canonical_amplitude_phase(amp, -x)
    end
    return (; amplitude, phase)
end
