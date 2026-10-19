const FACTORIAL = (
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0,
)

# Integer power by repeated multiplication: exact for the small exponents used here and free of
# the `0^0` and negative-base derivative pitfalls of `^`.
@inline function ipow(z, n::Integer)
    p = one(z)
    for _ in 1:n
        p *= z
    end
    return p
end

"""
    wigner_D_element(l, mp, m, ra, rb)

Element `D^l_{m′ m}` of the Wigner matrix of the rotor `R = w + x i + y j + z k`, written in the
spinor entries `ra = w + i z` and `rb = y + i x`. The closed-form sum has a fixed integer
bound set by `(l, m′, m)`: there is no square root, no case split on `|ra|` versus `|rb|` and no
branch on the rotor, so it is smooth at the identity and at `β = π`. Matches
`SphericalFunctions.D(R, l)[l][m′, m]`.
"""
function wigner_D_element(l::Integer, mp::Integer, m::Integer, ra::Complex, rb::Complex)
    0 <= l <= 5 || throw(ArgumentError("l must be in 0:5"))
    (abs(mp) <= l && abs(m) <= l) || throw(ArgumentError("need |m′|, |m| <= l"))
    ra_c = conj(ra)
    rb_c = conj(rb)
    acc = zero(ra)
    for s in max(0, m - mp):min(l + m, l - mp)
        den = FACTORIAL[l + m - s + 1] * FACTORIAL[s + 1] * FACTORIAL[mp - m + s + 1] * FACTORIAL[l - mp - s + 1]
        acc += ipow(ra_c, l + m - s) * ipow(ra, l - mp - s) * ipow(-rb, mp - m + s) * ipow(rb_c, s) / den
    end
    pref = sqrt(FACTORIAL[l + mp + 1] * FACTORIAL[l - mp + 1] * FACTORIAL[l + m + 1] * FACTORIAL[l - m + 1])
    return pref * acc
end

"""
    wigner_D(l, w, x, y, z)

The `(2l+1) × (2l+1)` Wigner matrix `D^l`, indexed `[m′ + l + 1, m + l + 1]`, of the unit
quaternion `(w, x, y, z)` (the caller is responsible for normalisation).
"""
function wigner_D(l::Integer, w, x, y, z)
    ra = complex(w, z)
    rb = complex(y, x)
    T = typeof(ra)
    D = Matrix{T}(undef, 2l + 1, 2l + 1)
    for m in (-l):l, mp in (-l):l
        D[mp + l + 1, m + l + 1] = wigner_D_element(l, mp, m, ra, rb)
    end
    return D
end
