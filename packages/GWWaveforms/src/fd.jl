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

"""
    hphc_fd(f::AbstractRange, params::NamedTuple, generator::Approximant) -> (; plus, cross)

Frequency-domain plus and cross polarisations on the uniform grid `f` (Hz), following bilby's
`lal_binary_black_hole` and the interface of [`AbstractGWWaveforms.hphc_fd`](@ref). `params` has bilby's
names and units: required `mass_1`, `mass_2` (detector-frame solar masses; other mass
parameters are converted by `GWParameterConversion.convert_parameters`),
`luminosity_distance` (Mpc), `theta_jn`, `phase`; optional (default 0) `a_1`, `a_2`, `tilt_1`,
`tilt_2`, `phi_12`, `phi_jl`, `eccentricity`.

`f` must start at 0 and have a positive step, as bilby assumes. The output is zero outside
`minimum_frequency ≤ f ≤ maximum_frequency`.

The parameter values are **not** checked: spins must be (anti-)aligned with the orbital angular
momentum, `eccentricity` zero, masses and distance positive, and so on. Run
[`validate_parameters`](@ref) first (once, on the host) for that.
`TaylorF2`, `IMRPhenomD` and `IMRPhenomXAS` support ForwardDiff and Mooncake.
`maximum_frequency` must be positive (`ArgumentError` otherwise).

The element type is `Complex{T}` with `T` the floating-point type promoted from the
parameters, so dual numbers (ForwardDiff) and Mooncake reverse mode pass through.
"""
function hphc_fd(f::AbstractRange, params::NamedTuple, g::Approximant)
    check_frequency_grid(f)
    T, p = promote_parameters(aligned_parameters(params))

    f_min, f_max = frequency_bounds(g, f)
    f_start = g.pn_amplitude_order != 0 ?
        f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min

    h = waveform_fd(g, T, f, p, f_start, f_min, f_max)

    # The models return h(f) for ι = 0; other orientations multiply by the inclination
    # dependence, hc ∝ -i hp (ChooseFDWaveform in LALSimInspiral).
    cfac = cos(p.inclination)
    pfac = 0.5 * (1.0 + cfac * cfac)
    re = real.(h)
    im = imag.(h)
    plus = complex.(pfac .* re, pfac .* im)
    cross = complex.(cfac .* im, .-(cfac .* re))
    return (; plus, cross)
end

# The inclination-independent kernels, selected by the generator type. They return zero outside
# `[f_min, f_max]`.
function waveform_fd(g::TaylorF2, T, f, p, f_start, f_min, f_max)
    lo, hi = waveform_slice(f, f_start, f_min, f_max)
    return taylorf2_fd_slice(T, length(f), lo, hi, Float64(step(f)), p, g)
end

function waveform_fd(g::IMRPhenomD, T, f, p, f_start, f_min, f_max)
    lo, hi = waveform_slice(f, f_start, f_min, f_max)
    return phenomd_fd_slice(T, length(f), lo, hi, Float64(step(f)), p, g, f_start, f_max)
end
