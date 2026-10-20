# Adapted from the MIT-licensed GWSurrogates Wigner-D polynomial (retained there for rotations).
# Precompute all coefficients at module load, never factorials or angle-dependent divisions
# during evaluation. Row l²+l+m-3 holds sqrt((2l+1)/4π) dˡ_{m,2} as a half-angle polynomial.
const HARMONIC_COEFFICIENTS = let
    a = zeros(32, 11)
    for l in 2:5, m in (-l):l
        pref = sqrt(Float64(factorial(l + m) * factorial(l - m)) * factorial(l + 2) * factorial(l - 2))
        for s in max(0, 2 - m):min(l + 2, l - m)
            n = m - 2 + 2s
            den = factorial(l + 2 - s) * factorial(s) * factorial(m - 2 + s) * factorial(l - m - s)
            a[l * l + l + m - 3, n + 1] = (-1)^(m - 2 + s) * pref / den * sqrt((2l + 1) / (4π))
        end
    end
    a
end

# Repeated multiplication avoids problematic derivatives of 0^0 and negative bases at poles.
function half_angle_powers(x)
    p = Vector{typeof(x)}(undef, 11)
    p[1] = one(x)
    for k in 2:11
        p[k] = p[k - 1] * x
    end
    return p
end

function harmonic_coefficient(l, m, cp, sp, phi)
    d = zero(eltype(cp))
    for n in 0:(2l)
        c = oftype(d, HARMONIC_COEFFICIENTS[l * l + l + m - 3, n + 1])
        d += c * cp[2l - n + 1] * sp[n + 1]
    end
    return d * cis(m * phi)
end

"""
    sylm_m2(l, m, theta, phi)

Spin-weight -2 spherical harmonic sqrt((2l+1)/4π) dˡ_{m,2}(theta) exp(im*m*phi),
for all signed modes with `2 ≤ l ≤ 5`. A real half-angle polynomial preserves finite
angle derivatives at both poles, without divisions by sin(theta).
"""
function sylm_m2(l::Integer, m::Integer, theta, phi)
    theta, phi = promote(float(theta), float(phi))
    return harmonic_coefficient(l, m, half_angle_powers(cos(theta / 2)), half_angle_powers(sin(theta / 2)), phi)
end

function mode_coefficients(labels, theta, phi)
    theta, phi = promote(float(theta), float(phi))
    cp = half_angle_powers(cos(theta / 2))
    sp = half_angle_powers(sin(theta / 2))
    return [harmonic_coefficient(l, m, cp, sp, phi) for (l, m) in labels]
end
