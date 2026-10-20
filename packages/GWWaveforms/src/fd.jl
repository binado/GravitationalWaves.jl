const MAX_FREQUENCY_ISCO_MESSAGE = "maximum_frequency must be positive: the default cutoff " *
    "is not exposed by hphc_fd; set an explicit upper frequency bound"

"""
    waveform_slice(f, f_start, f_min, f_max) -> (lo, hi)

Indices `lo:hi` (zero based, inclusive; `lo > hi` if empty) of the grid points where
the waveform is evaluated: `i ≥ ceil(f_start / Δf)`,
`i ≤ floor(f_max / Δf)` and `f_min ≤ f_i ≤ f_max`. Everything here comes from the grid and the
generator, not from the masses. Intrinsic model cutoffs are applied inside this interval.
"""
function waveform_slice(f::AbstractRange, f_start, f_min, f_max)
    f_max > 0 || throw(ArgumentError(MAX_FREQUENCY_ISCO_MESSAGE))
    deltaF = Float64(step(f))
    n = floor(Int, f_max / deltaF + 1)
    i_start = ceil(Int, f_start / deltaF)
    # Compare actual grid values to preserve exact frequency-bound semantics.
    first_in = 1
    while first_in <= length(f) && f[first_in] < f_min
        first_in += 1
    end
    last_in = length(f)
    while last_in >= 1 && f[last_in] > f_max
        last_in -= 1
    end
    return (max(i_start, first_in - 1), min(n, length(f), last_in) - 1)
end

# `XLALSimInspiralfLow2fStart`: the amplitude corrections of
# TaylorF2 start earlier than the requested frequency.
function f_low_to_f_start(f_low, amplitude_order::Integer)
    ampOrder = amplitude_order == -1 ? 6 : amplitude_order # MAX_NONPRECESSING_AMP_PN_ORDER
    return f_low * 2.0 / (ampOrder + 2)
end

# All parameters converted to one floating-point type, so the kernels are type stable and dual
# numbers in any parameter propagate through all of them.
function promote_parameters(p::NamedTuple)
    T = float(promote_type(map(typeof, values(p))...))
    return T, map(T, p)
end

# Convert a signed prefactor to the public nonnegative-amplitude, exp(+im * phase)
# convention without recovering phase from a complex strain or wrapping its analytic value.
function canonical_amplitude_phase(amplitude, phase)
    shifted = ifelse(amplitude < 0, phase + Float64(π), phase)
    return abs(amplitude), ifelse(iszero(amplitude), zero(shifted), shifted)
end

"""
    amplitude_phase(f::AbstractRange, params::NamedTuple, generator::Approximant)

Nonnegative amplitude and analytic unwrapped phase of the inclination-independent strain,
following the [`AbstractGWWaveforms.amplitude_phase`](@ref) convention `h = A * exp(im * Φ)`.
All three native models share the default quadrupole polarization assembly in `hphc_fd`.
Parameter values are trusted; run [`validate_parameters`](@ref) on the host first.
The grid and generator are fixed configuration, and `maximum_frequency` must be positive.
The real element type is the floating-point type promoted from the parameters, including duals.
"""
function amplitude_phase(f::AbstractRange, params::NamedTuple, g::Approximant)
    check_frequency_grid(f)
    T, p = promote_parameters(aligned_parameters(params))
    f_min, f_max = frequency_bounds(g, f)
    f_start = g.pn_amplitude_order != 0 ?
        f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min
    return amplitude_phase_fd(g, T, f, p, f_start, f_min, f_max)
end

function amplitude_phase_fd(g::TaylorF2, T, f, p, f_start, f_min, f_max)
    lo, hi = waveform_slice(f, f_start, f_min, f_max)
    return taylorf2_amplitude_phase_slice(T, length(f), lo, hi, Float64(step(f)), p, g)
end

function amplitude_phase_fd(g::IMRPhenomD, T, f, p, f_start, f_min, f_max)
    lo, hi = waveform_slice(f, f_start, f_min, f_max)
    return phenomd_amplitude_phase_slice(T, length(f), lo, hi, Float64(step(f)), p, g, f_start, f_max)
end
