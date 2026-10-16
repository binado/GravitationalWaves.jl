const REQUIRED_KEYS = (:mass_1, :mass_2, :luminosity_distance, :theta_jn, :phase)
const OPTIONAL_KEYS = (:a_1, :a_2, :tilt_1, :tilt_2, :phi_12, :phi_jl, :eccentricity)
const ALLOWED_KEYS = (REQUIRED_KEYS..., OPTIONAL_KEYS...)

"""
Binary parameters in bilby's names and units: detector-frame masses in solar masses,
luminosity distance in Mpc, angles in radians. Build from a `NamedTuple` via
`BBHParameters(nt)`, which rejects unknown keys.
"""
struct BBHParameters
    mass_1::Float64
    mass_2::Float64
    luminosity_distance::Float64
    theta_jn::Float64
    phase::Float64
    a_1::Float64
    a_2::Float64
    tilt_1::Float64
    tilt_2::Float64
    phi_12::Float64
    phi_jl::Float64
    eccentricity::Float64
end

BBHParameters(p::BBHParameters) = p

function BBHParameters(nt::NamedTuple)
    for k in keys(nt)
        k in ALLOWED_KEYS || throw(
            ArgumentError("unknown parameter `$k`; allowed keys are $(join(ALLOWED_KEYS, ", "))")
        )
    end
    for k in REQUIRED_KEYS
        haskey(nt, k) || throw(ArgumentError("missing required parameter `$k`"))
    end
    g(k) = Float64(get(nt, k, 0.0))
    return BBHParameters(
        g(:mass_1), g(:mass_2), g(:luminosity_distance), g(:theta_jn), g(:phase),
        g(:a_1), g(:a_2), g(:tilt_1), g(:tilt_2), g(:phi_12), g(:phi_jl), g(:eccentricity),
    )
end
