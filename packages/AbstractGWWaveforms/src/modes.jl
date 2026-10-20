"""
    HigherModesWaveformGenerator <: AbstractWaveformGenerator

Generator exposing fixed-frame spin-weight -2 modes through [`waveform_modes_td`](@ref)
or [`waveform_modes_fd`](@ref), inheriting polarization assembly in the supported domain.
No transform between domains is implied. Backends document their frame and reference phase.
"""
abstract type HigherModesWaveformGenerator <: AbstractWaveformGenerator end

"""
    waveform_modes_fd(f, params, generator) -> (; mode_indices, positive, negative)

Distance-normalized, inclination-independent complex modes on the nonnegative uniform grid
`f`. `mode_indices` is an ordered tuple of unique signed `(l,m)` labels with `2 ≤ l ≤ 5`.
Matrices have shape `(length(f), length(mode_indices))`; `positive[i,j]` is h̃ⱼ(f[i])
and `negative[i,j]` is h̃ⱼ(-f[i]), in the same row order. DC entries agree. Both branches
vanish outside configured frequency bounds and intrinsic model support. Reference phase
belongs in the modes by default; a backend may override [`mode_angles`](@ref).
"""
function waveform_modes_fd end

"""
    waveform_modes_td(generator, params; kwargs...) -> (; t, mode_indices, modes)

Fixed inertial/source-frame spherical modes in physical strain, excluding observer inclination.
`t` is a vector of seconds; complex `modes` has shape `(length(t), length(mode_indices))`.
Labels follow [`waveform_modes_fd`](@ref). The backend documents its reference-phase convention.
"""
function waveform_modes_td end

"""
    hphc_td(generator, params; kwargs...) -> (; t, hplus, hcross)

Time-domain polarizations in seconds and physical strain. Higher-mode generators evaluate
[`waveform_modes_td`](@ref) once and assemble H = Σ ₋₂Yₗₘ hₗₘ, hplus = real(H),
hcross = -imag(H), with angles from [`mode_angles`](@ref).
Mode labels and matrix dimensions are assumed to satisfy the [`waveform_modes_td`](@ref)
contract and are not revalidated during assembly.
"""
function hphc_td end

"""
    mode_angles(params, generator) -> (theta, phi)

Observer angles in the mode frame. Defaults to `(params.theta_jn, 0)`; the reference
phase is included in the modes. Backends may override this convention.
"""
mode_angles(params, ::HigherModesWaveformGenerator) = (params.theta_jn, zero(params.theta_jn))

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

"""
    check_mode_indices(mode_indices)

Validate a nonempty ordered tuple of unique integer `(l,m)` labels with `2 ≤ l ≤ 5`
and `|m| ≤ l`. Throws `ArgumentError` on invalid metadata.
"""
function check_mode_indices(labels)
    labels isa Tuple || throw(ArgumentError("mode_indices must be an ordered tuple"))
    isempty(labels) && throw(ArgumentError("mode_indices must be nonempty"))
    for (j, label) in enumerate(labels)
        label isa Tuple{Integer, Integer} || throw(ArgumentError("mode labels must be integer (l,m) pairs"))
        l, m = label
        2 <= l <= 5 && abs(m) <= l || throw(ArgumentError("need 2 <= l <= 5 and |m| <= l"))
        for k in 1:(j - 1)
            labels[k] == label && throw(ArgumentError("mode labels must be unique"))
        end
    end
    return nothing
end

function check_mode_matrix(labels, modes, n)
    modes isa AbstractMatrix{<:Complex} || throw(ArgumentError("modes must be a complex matrix"))
    size(modes) == (n, length(labels)) || throw(DimensionMismatch("mode matrix must have shape (sample, mode)"))
    Base.require_one_based_indexing(modes)
    return nothing
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
    check_mode_indices(((l, m),))
    theta, phi = promote(float(theta), float(phi))
    return harmonic_coefficient(l, m, half_angle_powers(cos(theta / 2)), half_angle_powers(sin(theta / 2)), phi)
end

"""
    mode_sum(modes, mode_indices, theta, phi)

Assemble complex strain H = Σ ₋₂Yₗₘ(theta,phi) hₗₘ from `(sample, mode)` complex
modes. Validate metadata once, compute angles once, and accumulate without weighted-array
intermediates. The azimuth `phi` is the actual mode-frame viewing angle.
"""
function mode_sum(modes::AbstractMatrix, labels, theta, phi)
    check_mode_indices(labels)
    check_mode_matrix(labels, modes, size(modes, 1))
    theta, phi = promote(float(theta), float(phi))
    coefficients = mode_coefficients(labels, theta, phi)
    return sum_modes(modes, coefficients)
end

function mode_coefficients(labels, theta, phi)
    theta, phi = promote(float(theta), float(phi))
    cp = half_angle_powers(cos(theta / 2))
    sp = half_angle_powers(sin(theta / 2))
    return [harmonic_coefficient(l, m, cp, sp, phi) for (l, m) in labels]
end

function sum_modes(modes, coefficients)
    T = promote_type(eltype(modes), eltype(coefficients))
    h = zeros(T, size(modes, 1))
    for j in eachindex(coefficients)
        y = coefficients[j]
        for k in eachindex(h)
            h[k] += y * modes[k, j]
        end
    end
    return h
end

function hphc_td(g::HigherModesWaveformGenerator, params::NamedTuple; kwargs...)
    r = waveform_modes_td(g, params; kwargs...)
    r.t isa AbstractVector{<:Real} || throw(ArgumentError("t must be a real vector of seconds"))
    theta, phi = mode_angles(params, g)
    h = sum_modes(r.modes, mode_coefficients(r.mode_indices, theta, phi))
    return (; t = r.t, hplus = real.(h), hcross = -imag.(h))
end

function hphc_fd(f::AbstractRange, params::NamedTuple, g::HigherModesWaveformGenerator)
    check_frequency_grid(f)
    r = waveform_modes_fd(f, params, g)
    check_mode_indices(r.mode_indices)
    check_mode_matrix(r.mode_indices, r.positive, length(f))
    check_mode_matrix(r.mode_indices, r.negative, length(f))
    theta, phi = mode_angles(params, g)
    coefficients = mode_coefficients(r.mode_indices, theta, phi)
    hp = sum_modes(r.positive, coefficients)
    hn = sum_modes(r.negative, coefficients)
    return (; plus = (hp .+ conj.(hn)) ./ 2, cross = (im / 2) .* (hp .- conj.(hn)))
end
