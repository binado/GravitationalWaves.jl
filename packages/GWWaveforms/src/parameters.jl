const REQUIRED_KEYS = (:mass_1, :mass_2, :luminosity_distance, :theta_jn, :phase)
const OPTIONAL_KEYS = (:a_1, :a_2, :tilt_1, :tilt_2, :phi_12, :phi_jl, :eccentricity)
# Any mass parameterisation is accepted; it is converted to mass_1 and mass_2.
const ALLOWED_KEYS = Tuple(unique((REQUIRED_KEYS..., OPTIONAL_KEYS..., MASS_KEYS...)))

"""
    aligned_parameters(params::NamedTuple) -> NamedTuple

Validate `params` (bilby's names and units, as for `LALWaveforms.hphc_fd`) and return
`(; mass_1, mass_2, luminosity_distance, inclination, phase, chi_1, chi_2)`: detector-frame
masses in solar masses, distance in Mpc, and the dimensionless spin components along the
orbital angular momentum. Element types are promoted from the inputs, so dual numbers pass through.

Unknown keys, a missing required key, spins that are not aligned or anti-aligned with the
orbital angular momentum, and a nonzero `eccentricity` throw `ArgumentError`: both models
are quasi-circular and aligned-spin, and LAL rejects such input too.
"""
function aligned_parameters(params::NamedTuple)
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
    spins = aligned_spins_to_lal_frame(nt.theta_jn, a_1, tilt_1, a_2, tilt_2)
    return (;
        mass_1 = nt.mass_1, mass_2 = nt.mass_2, luminosity_distance = nt.luminosity_distance,
        inclination = spins.iota, phase = nt.phase, chi_1 = spins.s1z, chi_2 = spins.s2z,
    )
end
