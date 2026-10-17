# Spin frame transform from bilby's parameterisation (angles relative to J and N) to the LAL
# source frame, as in `bilby_to_lalsimulation_spins` and
# `XLALSimInspiralTransformPrecessingNewInitialConditions` (LALSimInspiral.c, GPL-2.0-or-later).
#
# Written from the geometric definition, not transcribed: the sign and offset conventions
# (azimuth of L̂ about Ĵ is `phi_jl - π/2`, N̂ at azimuth `π/2 - phase` about L̂) and the
# 1PN |L| were fitted against the LAL library, and the LALWaveforms tests keep LAL as the
# oracle for this function.

"""
    is_aligned(a, tilt)

Whether a spin of magnitude `a` and angle `tilt` to the orbital angular momentum is aligned
or anti-aligned (or zero). Compares with the rounded double, not the `Irrational`, because
bilby tests `tilt in [0, np.pi]`.
"""
is_aligned(a::Number, tilt::Number) = (a == 0) | (tilt == 0) | (tilt == Float64(π))

"""
    aligned_spins_to_lal_frame(theta_jn, a_1, tilt_1, a_2, tilt_2) -> (; iota, s1x, s1y, s1z, s2x, s2y, s2z)

LAL-frame components for spins that are (anti-)aligned with L: the inclination is `theta_jn`
and only the z components are non-zero.
"""
function aligned_spins_to_lal_frame(
        theta_jn::Number, a_1::Number, tilt_1::Number, a_2::Number, tilt_2::Number
    )
    s1z = a_1 * cos(tilt_1)
    s2z = a_2 * cos(tilt_2)
    z = zero(s1z + s2z)
    return (; iota = theta_jn, s1x = z, s1y = z, s1z, s2x = z, s2y = z, s2z)
end

_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_norm(a) = sqrt(_dot(a, a))
_rotz(angle, v) = (cos(angle) * v[1] - sin(angle) * v[2], sin(angle) * v[1] + cos(angle) * v[2], v[3])
_roty(angle, v) = (cos(angle) * v[1] + sin(angle) * v[3], v[2], -sin(angle) * v[1] + cos(angle) * v[3])

"""
    transform_precessing_initial_conditions(theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2,
                                            mass_1, mass_2, f_ref, phase)
        -> (; iota, s1x, s1y, s1z, s2x, s2y, s2z)

Pure-Julia equivalent of `XLALSimInspiralTransformPrecessingNewInitialConditions`. Masses are
in solar masses and `f_ref` is the GW reference frequency in Hz; angles are in radians. The
returned spins are dimensionless components in the LAL source frame, where L̂ is the z axis
and the line of sight N̂ has inclination `iota` and azimuth `π/2 - phase`.

The orbital angular momentum is the 1PN expression at `v₀ = (π M f_ref)^(1/3)`. When N̂ is
(anti-)parallel to L̂ the azimuth reference is undefined; a fixed direction is used then.
"""
function transform_precessing_initial_conditions(
        theta_jn::Number, phi_jl::Number, tilt_1::Number, tilt_2::Number, phi_12::Number, a_1::Number,
        a_2::Number, mass_1::Number, mass_2::Number, f_ref::Number, phase::Number
    )
    # Promote values rather than constructing a scalar type: traced real numbers
    # need not subtype Real or support type-level floating-point construction.
    exemplar = float(
        first(
            promote(
                theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2, mass_1, mass_2, f_ref, phase,
            )
        )
    )
    zero_T = zero(exemplar)
    one_T = one(exemplar)
    pi_T = oftype(exemplar, π)

    # Starting frame: L̂ along z, spin 1 in the x-z plane, spin 2 at azimuth phi_12 from it.
    s1hat = (sin(tilt_1), zero_T, cos(tilt_1))
    s2hat = (sin(tilt_2) * cos(phi_12), sin(tilt_2) * sin(phi_12), cos(tilt_2))

    total = mass_1 + mass_2
    eta = mass_1 * mass_2 / (total * total)
    v0 = cbrt(total * MTSUN_SI * π * f_ref)
    l_mag = total * total * eta / v0 * (1 + v0 * v0 * (oftype(exemplar, 3) / 2 + eta / 6))

    j = (
        mass_1^2 * a_1 * s1hat[1] + mass_2^2 * a_2 * s2hat[1],
        mass_1^2 * a_1 * s1hat[2] + mass_2^2 * a_2 * s2hat[2],
        l_mag + mass_1^2 * a_1 * s1hat[3] + mass_2^2 * a_2 * s2hat[3],
    )
    j_norm = _norm(j)
    theta_0 = acos(j[3] / j_norm)
    phi_0 = atan(j[2], j[1])

    # J frame: Ĵ along z and N̂ in the x-z plane; L̂ sits at azimuth `phi_jl - π/2` on its cone.
    to_j(v) = _rotz(phi_jl - pi_T / 2 - pi_T, _roty(-theta_0, _rotz(-phi_0, v)))
    l_hat = to_j((zero_T, zero_T, one_T))
    s1_j = to_j(s1hat)
    s2_j = to_j(s2hat)
    n_hat = (sin(theta_jn), zero_T, cos(theta_jn))

    iota = atan(_norm(_cross(n_hat, l_hat)), _dot(n_hat, l_hat))

    # LAL source frame: z = L̂, with N̂'s projection at azimuth π/2 - phase.
    n_perp = n_hat .- _dot(n_hat, l_hat) .* l_hat
    n_perp_sq = _dot(n_perp, n_perp)
    use_projection = sqrt(n_perp_sq) > 4 * eps(one_T)
    # Select safe radicands before sqrt/division. The unused projection can be
    # zero, and the unused fallback can be singular when L is parallel to y.
    n_perp_norm = sqrt(ifelse(use_projection, n_perp_sq, one_T))
    y_cross = _cross((zero_T, one_T, zero_T), l_hat)
    y_cross_norm = sqrt(ifelse(use_projection, one_T, _dot(y_cross, y_cross)))
    n_perp_hat = map(n_perp, y_cross) do n, y
        ifelse(use_projection, n / n_perp_norm, y / y_cross_norm)
    end
    azimuth = pi_T / 2 - phase
    x_hat = cos(azimuth) .* n_perp_hat .- sin(azimuth) .* _cross(l_hat, n_perp_hat)
    y_hat = _cross(l_hat, x_hat)

    components(s, a) = (a * _dot(s, x_hat), a * _dot(s, y_hat), a * _dot(s, l_hat))
    s1 = components(s1_j, a_1)
    s2 = components(s2_j, a_2)
    return (; iota, s1x = s1[1], s1y = s1[2], s1z = s1[3], s2x = s2[1], s2y = s2[2], s2z = s2[3])
end

"""
    spins_to_lal_frame(p::NamedTuple, f_ref) -> (; iota, s1x, s1y, s1z, s2x, s2y, s2z)

bilby's `bilby_to_lalsimulation_spins` for a parameter `NamedTuple` with detector-frame
`mass_1`, `mass_2` (solar masses), `theta_jn`, `phase`, and optional (default 0) `a_1`, `a_2`,
`tilt_1`, `tilt_2`, `phi_12`, `phi_jl`. Aligned spins select the exact aligned result; otherwise
[`transform_precessing_initial_conditions`](@ref) is used.
"""
function spins_to_lal_frame(p::NamedTuple, f_ref::Number)
    a_1, a_2 = get(p, :a_1, 0.0), get(p, :a_2, 0.0)
    tilt_1, tilt_2 = get(p, :tilt_1, 0.0), get(p, :tilt_2, 0.0)
    aligned = is_aligned(a_1, tilt_1) & is_aligned(a_2, tilt_2)
    aligned_result = aligned_spins_to_lal_frame(p.theta_jn, a_1, tilt_1, a_2, tilt_2)
    # Eager numeric selection evaluates both paths. Substitute the nonsingular
    # precessing regression configuration before evaluating the inactive path.
    inputs = (
        p.theta_jn, get(p, :phi_jl, 0.0), tilt_1, tilt_2, get(p, :phi_12, 0.0), a_1, a_2,
        p.mass_1, p.mass_2, f_ref, p.phase,
    )
    safe = (0.4, 2.5, 0.7, 1.9, 1.1, 0.6, 0.4, 36.0, 29.0, 50.0, 1.3)
    selected = map((x, fallback) -> ifelse(aligned, oftype(float(x), fallback), x), inputs, safe)
    precessing_result = transform_precessing_initial_conditions(selected...)
    return map((a, b) -> ifelse(aligned, a, b), aligned_result, precessing_result)
end
