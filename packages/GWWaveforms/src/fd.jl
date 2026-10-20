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
