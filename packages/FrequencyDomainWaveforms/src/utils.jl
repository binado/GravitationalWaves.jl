const MAX_FREQUENCY_ISCO_MESSAGE = "maximum_frequency must be positive: the default cutoff " *
    "is not exposed by hphc_fd; set an explicit upper frequency bound"

# `XLALSimInspiralfLow2fStart`: the amplitude corrections of
# TaylorF2 start earlier than the requested frequency.
function f_low_to_f_start(f_low, amplitude_order::Integer)
    amp_order = amplitude_order == -1 ? 6 : amplitude_order # MAX_NONPRECESSING_AMP_PN_ORDER
    return f_low * 2.0 / (amp_order + 2)
end

# Lowest frequency at which the waveform is evaluated, given the generator's amplitude order.
start_frequency(g, f_min) = g.pn_amplitude_order != 0 ? f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min

# All parameters converted to one floating-point type, so the kernels are type stable and dual
# numbers in any parameter propagate through all of them.
function promote_parameters(params::NamedTuple)
    T = float(promote_type(map(typeof, values(params))...))
    return T, map(T, params)
end

# Convert a signed prefactor to the public nonnegative-amplitude, exp(+im * phase)
# convention without recovering phase from a complex strain or wrapping its analytic value.
function canonical_amplitude_phase(amplitude, phase)
    shifted = ifelse(amplitude < 0, phase + Float64(π), phase)
    return abs(amplitude), ifelse(iszero(amplitude), zero(shifted), shifted)
end

# One-based indices into the requested grid; no relationship between index and frequency
# is assumed. Grid and generator validation runs separately on the host.
function frequency_slice(f::AbstractRange{Float64}, f_start, f_min, f_max)
    return (searchsortedfirst(f, max(f_start, f_min)), searchsortedlast(f, f_max))
end
function check_pn_orders(; pn_spin_order = 0, pn_phase_order = 0, pn_amplitude_order = 0)
    -1 <= pn_spin_order <= 7 || throw(ArgumentError("invalid pn_spin_order $pn_spin_order"))
    -1 <= pn_phase_order <= 7 || throw(ArgumentError("invalid pn_phase_order $pn_phase_order"))
    # LAL has no 1PN-relative amplitude correction.
    (pn_amplitude_order == -1 || pn_amplitude_order == 0 || 2 <= pn_amplitude_order <= 7) ||
        throw(ArgumentError("invalid pn_amplitude_order $pn_amplitude_order"))
    return nothing
end

maximum_or_nothing(f::Nothing) = nothing
maximum_or_nothing(f::Real) = Float64(f)

# LAL uses the starting frequency when no reference frequency is specified.
phenom_reference_frequency(g, f_start) = g.reference_frequency == 0.0 ? f_start : g.reference_frequency

# Shared strain scale; keep LAL's multiplication order for numerical fidelity.
function phenom_amplitude_scale(m1, m2, luminosity_distance)
    M = m1 + m2
    distance = luminosity_distance * 1.0e6 * PC_SI
    return 2.0 * sqrt(5.0 / (64.0 * Float64(π))) * M * MRSUN_SI * M * MTSUN_SI / distance
end

# Host-only physical checks; remnant and frequency restrictions belong to each model.
function check_phenom_parameters(params::NamedTuple)
    (params.mass_1 > 0 && params.mass_2 > 0) || throw(ArgumentError("masses must be positive"))
    params.luminosity_distance > 0 || throw(ArgumentError("luminosity distance must be positive"))
    (abs(params.chi_1) <= 1 && abs(params.chi_2) <= 1) ||
        throw(ArgumentError("spins outside the range [-1, 1] are not supported"))
    return nothing
end
