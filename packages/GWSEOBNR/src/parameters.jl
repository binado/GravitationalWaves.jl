primal(x::Real) = x
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))

const MASS_KEYS = GWParameterConversion.MASS_KEYS
const REQUIRED_KEYS = (:mass_1, :mass_2, :luminosity_distance, :theta_jn, :phase)
const OPTIONAL_KEYS = (:a_1, :a_2, :tilt_1, :tilt_2, :phi_12, :phi_jl, :eccentricity)
const ALLOWED_KEYS = Tuple(unique((REQUIRED_KEYS..., OPTIONAL_KEYS..., MASS_KEYS...)))

function aligned_parameters(params::NamedTuple)
    nt = convert_parameters(params, (:mass_1, :mass_2); check_consistency = false)
    return (;
        nt.mass_1, nt.mass_2, nt.luminosity_distance,
        inclination = nt.theta_jn, nt.phase,
        chi_1 = get(nt, :a_1, 0.0) * cos(get(nt, :tilt_1, 0.0)),
        chi_2 = get(nt, :a_2, 0.0) * cos(get(nt, :tilt_2, 0.0)),
    )
end

function validate_parameters(f::AbstractRange, params::NamedTuple, g::SEOBNRv5HM)
    check_frequency_grid(f)
    for k in keys(params)
        k in ALLOWED_KEYS || throw(ArgumentError("unknown parameter `$k`"))
        isfinite(getproperty(params, k)) || throw(ArgumentError("parameter `$k` must be finite"))
    end
    nt = convert_parameters(params, (:mass_1, :mass_2))
    for k in REQUIRED_KEYS
        haskey(nt, k) || throw(ArgumentError("missing required parameter `$k`"))
    end
    for i in (1, 2)
        a = get(nt, Symbol(:a_, i), 0.0)
        tilt = get(nt, Symbol(:tilt_, i), 0.0)
        0 <= a <= 0.6 || throw(ArgumentError("spin magnitudes must lie in [0, 0.6]"))
        is_aligned(primal(a), Float64(primal(tilt))) || throw(ArgumentError("only aligned or anti-aligned spins are supported"))
    end
    iszero(get(nt, :eccentricity, 0.0)) || throw(ArgumentError("eccentricity must be zero"))
    p = aligned_parameters(nt)
    40 <= p.mass_1 + p.mass_2 <= 100 || throw(ArgumentError("total detector-frame mass must lie in [40, 100] solar masses"))
    p.mass_2 > 0 && 1 <= p.mass_1 / p.mass_2 <= 4 || throw(ArgumentError("ordered mass ratio mass_1/mass_2 must lie in [1, 4]"))
    p.luminosity_distance > 0 || throw(ArgumentError("luminosity_distance must be positive"))
    0 <= p.inclination <= π || throw(ArgumentError("theta_jn must lie in [0, π]"))
    g.minimum_frequency < last(f) || throw(ArgumentError("grid must extend above minimum_frequency"))
    isnothing(g.maximum_frequency) || g.maximum_frequency <= last(f) ||
        throw(ArgumentError("maximum_frequency exceeds the grid"))
    return nothing
end
