using BackgroundCosmology: AbstractCosmology, FlatLCDM, lookback_time, age, E, H0, Ωm

# Hubble time 1/H₀ in Gyr. BackgroundCosmology's own version is internal, so it is restated
# here; the constants match its conversions (Mpc = 3.085677581e22 m, Julian-year Gyr).
const _METERS_PER_MPC = 3.085677581e22
const _SECONDS_PER_GYR = JULIAN_YEAR_SEC * 1.0e9

_hubble_time_gyr(H0_km_s_mpc::Real) = _METERS_PER_MPC / (1000 * H0_km_s_mpc) / _SECONDS_PER_GYR

"""
    _redshift_at_lookback_time(c, t) -> Real

Redshift at lookback time `t` (Gyr), the inverse of `lookback_time(c, ·)`. Returns `Inf`
at and beyond the age of the universe.

For [`FlatLCDM`] this is the closed form

    1 + z = [√(ΩΛ/Ωm) / sinh A]^(2/3),   A = asinh√(ΩΛ/Ωm) - (3√ΩΛ / 2) t / t_H.

Other models use a fixed number of Newton iterations (so `ForwardDiff` duals propagate) in
`w = (1+z)^(-3/2)`, in which the cosmic age is nearly linear.
"""
function _redshift_at_lookback_time(c::FlatLCDM, t::Real)
    sqrt_ΩΛ = sqrt(1 - Ωm(c))
    ratio = sqrt_ΩΛ / sqrt(Ωm(c))
    angle = asinh(ratio) - 1.5 * sqrt_ΩΛ * t / _hubble_time_gyr(H0(c))
    angle > 0 || return oftype(angle, Inf)
    return (ratio / sinh(angle))^(2 / 3) - 1
end

const _NEWTON_ITERATIONS = 12

function _redshift_at_lookback_time(c::AbstractCosmology, t::Real)
    age0 = age(c, 0)
    target = age0 - t
    target > 0 || return oftype(target, Inf)
    t_H = _hubble_time_gyr(H0(c))
    # age(w) is concave and increasing in w, nearly linear for small w; start from the
    # linear guess through the origin and keep w inside (0, 1].
    w = clamp(target / age0, 1.0e-12, 1.0) * one(target)
    for _ in 1:_NEWTON_ITERATIONS
        z = w^(-2 / 3) - 1
        residual = age(c, z) - target
        slope = (2 / 3) * t_H / (w * E(c, z))
        w = clamp(w - residual / slope, 1.0e-12, 1.0)
    end
    return w^(-2 / 3) - 1
end

"""
    _log_gw_em_ratio(z, xi_0, xi_n)

Log of the modified-propagation GW-to-EM luminosity-distance ratio
`xi_0 + (1 - xi_0) (1 + z)^(-xi_n)`; identically 0 for `xi_0 = 1`.
"""
function _log_gw_em_ratio(z::Real, xi_0::Real, xi_n::Real)
    return log(xi_0 + (1 - xi_0) * exp(-(xi_n * log1p(z))))
end
