const FACTORIAL = (
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0,
)

# Powers `z^0, …, z^n` of the spinor entries, stored at index `k + 1`. Every term of the closed-form
# sum needs integer powers of `ra`, `conj(ra)`, `-rb` and `conj(rb)` (total degree `2l`); building them
# once per rotor instead of once per `(l, m′, m, s)` term removes most of the complex arithmetic,
# which is the expensive part when the entries carry dual numbers. Integer powers by repeated
# multiplication are exact and free of the `0^0` and negative-base derivative pitfalls of `^`.
struct SpinorPowers{T}
    ra::Vector{T}
    ra_c::Vector{T}
    nrb::Vector{T}
    rb_c::Vector{T}
end

SpinorPowers(::Type{T}, degree::Integer) where {T} = SpinorPowers{T}(
    Vector{T}(undef, degree + 1), Vector{T}(undef, degree + 1), Vector{T}(undef, degree + 1), Vector{T}(undef, degree + 1),
)

function fill_powers!(p::SpinorPowers, ra, rb)
    p.ra[1] = one(ra)
    p.ra_c[1] = one(ra)
    p.nrb[1] = one(rb)
    p.rb_c[1] = one(rb)
    nrb = -rb
    rb_c = conj(rb)
    ra_c = conj(ra)
    for k in 2:length(p.ra)
        p.ra[k] = p.ra[k - 1] * ra
        p.ra_c[k] = p.ra_c[k - 1] * ra_c
        p.nrb[k] = p.nrb[k - 1] * nrb
        p.rb_c[k] = p.rb_c[k - 1] * rb_c
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
    p = fill_powers!(SpinorPowers(typeof(ra), 2l), ra, rb)
    return wigner_D_element(l, mp, m, p)
end

# Same element from precomputed `SpinorPowers` of degree at least `2l`; `l`, `m′`, `m` are not validated.
function wigner_D_element(l::Integer, mp::Integer, m::Integer, p::SpinorPowers)
    acc = zero(eltype(p.ra))
    for s in max(0, m - mp):min(l + m, l - mp)
        den = FACTORIAL[l + m - s + 1] * FACTORIAL[s + 1] * FACTORIAL[mp - m + s + 1] * FACTORIAL[l - mp - s + 1]
        acc += p.ra_c[l + m - s + 1] * p.ra[l - mp - s + 1] * p.nrb[mp - m + s + 1] * p.rb_c[s + 1] / den
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
    0 <= l <= 5 || throw(ArgumentError("l must be in 0:5"))
    ra = complex(w, z)
    rb = complex(y, x)
    T = typeof(ra)
    p = fill_powers!(SpinorPowers(T, 2l), ra, rb)
    D = Matrix{T}(undef, 2l + 1, 2l + 1)
    for m in (-l):l, mp in (-l):l
        D[mp + l + 1, m + l + 1] = wigner_D_element(l, mp, m, p)
    end
    return D
end
