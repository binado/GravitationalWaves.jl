using Distributions
using Random

"""
    UniformCosine()

Polar angle `θ = acos(U)`, `U ~ Uniform(-1, 1)`, of a direction uniform on the sphere (the
law of the inclination of an isotropic population). Density `sin(θ) / 2` on `[0, π]`,
CDF `(1 - cos θ) / 2` (saturating outside), quantile `acos(1 - 2q)`.
"""
struct UniformCosine <: ContinuousUnivariateDistribution end

Base.minimum(::UniformCosine) = 0.0
Base.maximum(::UniformCosine) = float(π)
Distributions.insupport(d::UniformCosine, v::Real) = minimum(d) <= v <= maximum(d)
Base.eltype(::Type{UniformCosine}) = Float64

function Distributions.logpdf(d::UniformCosine, v::Real)
    insupport(d, v) || return oftype(float(v), -Inf)
    return log(sin(v)) - log(2)
end

Distributions.cdf(::UniformCosine, v::Real) = (1 - cos(clamp(v, 0, π))) / 2

Distributions.quantile(::UniformCosine, q::Real) = acos(1 - 2q)

Random.rand(rng::AbstractRNG, ::UniformCosine) = acos(2 * rand(rng) - 1)
