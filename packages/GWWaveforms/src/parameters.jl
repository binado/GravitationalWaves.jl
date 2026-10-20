const REQUIRED_KEYS = (:mass_1, :mass_2, :luminosity_distance, :theta_jn, :phase)
const OPTIONAL_KEYS = (:a_1, :a_2, :tilt_1, :tilt_2, :phi_12, :phi_jl, :eccentricity)
# Any mass parameterisation is accepted; it is converted to mass_1 and mass_2.
const ALLOWED_KEYS = Tuple(unique((REQUIRED_KEYS..., OPTIONAL_KEYS..., MASS_KEYS...)))

"""
    aligned_parameters(params::NamedTuple) -> NamedTuple

Convert `params` (bilby's names and units, as for `LALWaveforms.hphc_fd`) to
`(; mass_1, mass_2, luminosity_distance, inclination, phase, chi_1, chi_2)`: detector-frame
masses in solar masses, distance in Mpc, and the dimensionless spin components along the
orbital angular momentum. Element types are promoted from the inputs, so dual numbers pass through.

The values are not checked; see
[`validate_parameters`](@ref) for that. The spin components are `a cos(tilt)`, which is exact
for the (anti-)aligned spins the models support.
"""
function aligned_parameters(params::NamedTuple)
    nt = convert_parameters(params, (:mass_1, :mass_2); check_consistency = false)
    spins = aligned_spins_to_lal_frame(
        nt.theta_jn, get(nt, :a_1, 0.0), get(nt, :tilt_1, 0.0), get(nt, :a_2, 0.0), get(nt, :tilt_2, 0.0)
    )
    return (;
        mass_1 = nt.mass_1, mass_2 = nt.mass_2, luminosity_distance = nt.luminosity_distance,
        inclination = spins.iota, phase = nt.phase, chi_1 = spins.s1z, chi_2 = spins.s2z,
    )
end

# The checks on `params` alone: names, redundant mass parameters, and the model assumptions.
function check_parameters(params::NamedTuple)
    for k in keys(params)
        k in ALLOWED_KEYS || throw(
            ArgumentError("unknown parameter `$k`; allowed keys are $(join(ALLOWED_KEYS, ", "))")
        )
    end
    nt = convert_parameters(params, (:mass_1, :mass_2))
    for k in REQUIRED_KEYS
        haskey(nt, k) || throw(ArgumentError("missing required parameter `$k`"))
    end
    a_1, a_2 = get(nt, :a_1, 0.0), get(nt, :a_2, 0.0)
    tilt_1, tilt_2 = get(nt, :tilt_1, 0.0), get(nt, :tilt_2, 0.0)
    (is_aligned(a_1, tilt_1) && is_aligned(a_2, tilt_2)) || throw(
        ArgumentError(
            "only aligned or anti-aligned spins are supported (tilt_1 and tilt_2 must be 0 or π, " *
                "or the spin magnitude zero)"
        )
    )
    iszero(get(nt, :eccentricity, 0.0)) ||
        throw(ArgumentError("eccentric orbits are not supported (eccentricity must be 0)"))
    return nothing
end

"""
    validate_parameters(f::AbstractRange, params::NamedTuple, generator::Approximant) -> nothing

Everything [`hphc_fd`](@ref) assumes but does not check: the grid, the parameter names, redundant mass parameters, spins that are not aligned
or anti-aligned with the orbital angular momentum, a nonzero `eccentricity`, and the ranges of the
masses, distance and frequencies. Throws `ArgumentError` (`InconsistentParametersError` for
redundant mass parameters that disagree) and returns `nothing` otherwise. All models are
quasi-circular and aligned-spin, and LAL rejects such input too.

Call it on the host with the parameters of one draw, before differentiating or sampling `hphc_fd`.
"""
function validate_parameters(f::AbstractRange, params::NamedTuple, g::Approximant)
    check_frequency_grid(f)
    check_parameters(params)
    p = aligned_parameters(params)
    f_min, f_max = frequency_bounds(g, f)
    f_start = g.pn_amplitude_order != 0 ?
        f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min
    check_waveform_inputs(g, step(f), p, f_start, f_max)
    return nothing
end
