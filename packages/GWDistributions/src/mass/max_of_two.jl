using Distributions
using Random

"""
    MaxOfTwoNormals(loc = 0, scale = 1)

Law of `max(X, Y)` for i.i.d. `X, Y ~ Normal(loc, scale)`: density `2 φ(x) Φ(x)`, CDF
`Φ(x)²`, quantile `μ + σ Φ⁻¹(√q)`. Paired with a `TruncatedNormal(…, upper = m₁)`
secondary it gives the ordered joint `2 N(m₁) N(m₂)` on `m₁ ≥ m₂`.
"""
struct MaxOfTwoNormals{T <: Real} <: ContinuousUnivariateDistribution
    loc::T
    scale::T
    function MaxOfTwoNormals{T}(loc::T, scale::T) where {T <: Real}
        scale > 0 || throw(ArgumentError("scale must be positive"))
        return new{T}(loc, scale)
    end
end

function MaxOfTwoNormals(loc::Real = 0, scale::Real = 1)
    T = promote_type(typeof(float(loc)), typeof(float(scale)))
    return MaxOfTwoNormals{T}(T(loc), T(scale))
end

Base.minimum(::MaxOfTwoNormals{T}) where {T} = T(-Inf)
Base.maximum(::MaxOfTwoNormals{T}) where {T} = T(Inf)
Distributions.insupport(::MaxOfTwoNormals, v::Real) = isfinite(v)
Base.eltype(::Type{<:MaxOfTwoNormals{T}}) where {T} = T

_normal(d::MaxOfTwoNormals) = Normal(d.loc, d.scale)

function Distributions.logpdf(d::MaxOfTwoNormals, v::Real)
    n = _normal(d)
    return log(2) + logpdf(n, v) + logcdf(n, v)
end

Distributions.cdf(d::MaxOfTwoNormals, v::Real) = cdf(_normal(d), v)^2

Distributions.quantile(d::MaxOfTwoNormals, q::Real) = quantile(_normal(d), sqrt(q))

function Random.rand(rng::AbstractRNG, d::MaxOfTwoNormals)
    return max(rand(rng, _normal(d)), rand(rng, _normal(d)))
end

"""
    MaxOfTwoUniforms(low = 0, width = 1)

Law of `max(X, Y)` for i.i.d. `X, Y ~ Uniform(low, low + width)`, i.e. a scaled
`Beta(2, 1)`: density `2 (x - low) / width²`, CDF `((x - low) / width)²`, quantile
`low + width √q`. Paired with a `Uniform(low, m₁)` secondary it gives the ordered joint
of two i.i.d. uniforms.
"""
struct MaxOfTwoUniforms{T <: Real} <: ContinuousUnivariateDistribution
    low::T
    width::T
    function MaxOfTwoUniforms{T}(low::T, width::T) where {T <: Real}
        width > 0 || throw(ArgumentError("width must be positive"))
        return new{T}(low, width)
    end
end

function MaxOfTwoUniforms(low::Real = 0, width::Real = 1)
    T = promote_type(typeof(float(low)), typeof(float(width)))
    return MaxOfTwoUniforms{T}(T(low), T(width))
end

Base.minimum(d::MaxOfTwoUniforms) = d.low
Base.maximum(d::MaxOfTwoUniforms) = d.low + d.width
Distributions.insupport(d::MaxOfTwoUniforms, v::Real) = minimum(d) <= v <= maximum(d)
Base.eltype(::Type{<:MaxOfTwoUniforms{T}}) where {T} = T

function Distributions.logpdf(d::MaxOfTwoUniforms, v::Real)
    insupport(d, v) || return oftype(float(v), -Inf)
    return log(2) + log(v - d.low) - 2 * log(d.width)
end

function Distributions.cdf(d::MaxOfTwoUniforms, v::Real)
    return clamp((v - d.low) / d.width, 0, 1)^2
end

Distributions.quantile(d::MaxOfTwoUniforms, q::Real) = d.low + d.width * sqrt(q)

Random.rand(rng::AbstractRNG, d::MaxOfTwoUniforms) = quantile(d, rand(rng))
