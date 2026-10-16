using Distributions
using Random

# Below this magnitude the closed forms of `expm1(x) / x` and `log1p(y) / y` cancel; the
# truncation error of the series is O(x³) ≈ 1e-19 there.
const _SERIES_THRESHOLD = 1.0e-6

"""`expm1(x) / x`, equal to 1 at 0 and smooth through it."""
@inline function _expm1_ratio(x::Real)
    abs(x) < _SERIES_THRESHOLD && return 1 + x / 2 + x^2 / 6
    return expm1(x) / x
end

"""`log1p(y) / y`, equal to 1 at 0 and smooth through it."""
@inline function _log1p_ratio(y::Real)
    abs(y) < _SERIES_THRESHOLD && return 1 - y / 2 + y^2 / 3
    return log1p(y) / y
end

"""
    PowerLawDelay(alpha, low, high)

Delay-time distribution `p(τ) ∝ τ^alpha` on `[low, high]` (`low > 0`), written in
`β = 1 + alpha` through `expm1(x)/x` and `log1p(y)/y` so that `cdf`, `quantile` and their
derivatives in `alpha` stay finite and smooth through the canonical `alpha = -1`. With
`ℓ = log(high/low)` and `s = log(τ/low)`:

    F(τ)  = s h(βs) / (ℓ h(βℓ)),        F⁻¹(u) = low exp[g(y) u ℓ h(βℓ)],  y = u β ℓ h(βℓ)

for `h(x) = expm1(x)/x` and `g(y) = log1p(y)/y`.
"""
struct PowerLawDelay{T <: Real} <: ContinuousUnivariateDistribution
    alpha::T
    low::T
    high::T
    function PowerLawDelay{T}(alpha::T, low::T, high::T) where {T <: Real}
        0 < low < high || throw(ArgumentError("PowerLawDelay requires 0 < low < high"))
        return new{T}(alpha, low, high)
    end
end

function PowerLawDelay(alpha::Real, low::Real, high::Real)
    T = promote_type(typeof(float(alpha)), typeof(float(low)), typeof(float(high)))
    return PowerLawDelay{T}(T(alpha), T(low), T(high))
end

Base.minimum(d::PowerLawDelay) = d.low
Base.maximum(d::PowerLawDelay) = d.high
Distributions.insupport(d::PowerLawDelay, v::Real) = d.low <= v <= d.high
Base.eltype(::Type{<:PowerLawDelay{T}}) where {T} = T

@inline _delay_shape(d::PowerLawDelay) = (1 + d.alpha, log(d.high / d.low))

function Distributions.logpdf(d::PowerLawDelay, v::Real)
    insupport(d, v) || return oftype(float(v) * one(d.alpha), -Inf)
    β, ℓ = _delay_shape(d)
    log_norm = β * log(d.low) + log(ℓ) + log(_expm1_ratio(β * ℓ))
    return d.alpha * log(v) - log_norm
end

function Distributions.cdf(d::PowerLawDelay, v::Real)
    β, ℓ = _delay_shape(d)
    s = log(max(v, d.low) / d.low)
    return clamp((s * _expm1_ratio(β * s)) / (ℓ * _expm1_ratio(β * ℓ)), 0, 1)
end

function Distributions.quantile(d::PowerLawDelay, q::Real)
    β, ℓ = _delay_shape(d)
    scaled = ℓ * _expm1_ratio(β * ℓ)
    return d.low * exp(_log1p_ratio(q * β * scaled) * q * scaled)
end

Random.rand(rng::AbstractRNG, d::PowerLawDelay) = quantile(d, rand(rng))
