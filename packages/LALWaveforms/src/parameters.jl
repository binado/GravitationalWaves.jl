const REQUIRED_KEYS = (:mass_1, :mass_2, :luminosity_distance, :theta_jn, :phase)
const OPTIONAL_KEYS = (:a_1, :a_2, :tilt_1, :tilt_2, :phi_12, :phi_jl, :eccentricity)
# Any mass parameterisation is accepted; `bilby_parameters` derives mass_1 and mass_2 from it.
const ALLOWED_KEYS = Tuple(unique((REQUIRED_KEYS..., OPTIONAL_KEYS..., MASS_KEYS...)))

"""
    bilby_parameters(params::NamedTuple) -> NamedTuple

Validate `params` and return the full `Float64` `NamedTuple` that the LAL calls use, with
bilby's names and units: detector-frame masses in solar masses, luminosity distance in Mpc,
angles in radians. Optional parameters default to 0. Unknown keys throw `ArgumentError`.

`mass_1` and `mass_2` may be given through any other mass parameters that
`GWParameterConversion.convert_parameters` can convert (`chirp_mass` and `mass_ratio`, ...);
redundant masses must agree, else `InconsistentParametersError` is thrown.
"""
function bilby_parameters(params::NamedTuple)
    for k in keys(params)
        k in ALLOWED_KEYS || throw(
            ArgumentError("unknown parameter `$k`; allowed keys are $(join(ALLOWED_KEYS, ", "))")
        )
    end
    nt = convert_parameters(params, (:mass_1, :mass_2))
    for k in REQUIRED_KEYS
        haskey(nt, k) || throw(ArgumentError("missing required parameter `$k`"))
    end
    g(k) = Float64(get(nt, k, 0.0))
    return (;
        mass_1 = g(:mass_1), mass_2 = g(:mass_2), luminosity_distance = g(:luminosity_distance),
        theta_jn = g(:theta_jn), phase = g(:phase), a_1 = g(:a_1), a_2 = g(:a_2),
        tilt_1 = g(:tilt_1), tilt_2 = g(:tilt_2), phi_12 = g(:phi_12), phi_jl = g(:phi_jl),
        eccentricity = g(:eccentricity),
    )
end
