using Distributions
using Random

"""
    InterpolatedDistribution(x, y)

Tabulated unnormalized density `y` on the increasing grid `x`, with a linear interpolant.
The density is `interp(y) / norm`, where [`normalizer`](@ref) `= ∫ y dx` by the trapezoid
rule, and the CDF is the trapezoid cumulative integral at the nodes (linear in between),
so [`logpdf`](@ref), `cdf`, `quantile` and `rand` all derive from the same table. The
support is the closed interval `[x[1], x[end]]`.

Operation order follows `astrogwb`: `logpdf = log(interp(y) / norm)`, which gives `-Inf`
off the grid and wherever `y = 0`.
"""
struct InterpolatedDistribution{
        TX <: AbstractVector, TY <: AbstractVector, TN <: Real, TC <: AbstractVector,
    } <: ContinuousUnivariateDistribution
    x::TX
    y::TY
    norm::TN
    cdf_grid::TC
end

function InterpolatedDistribution(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(
        DimensionMismatch("InterpolatedDistribution grid and density must have the same length")
    )
    length(x) >= 2 || throw(
        ArgumentError("InterpolatedDistribution requires at least two grid points")
    )
    norm = _trapz(y, x)
    cdf = _cumtrapz(y ./ norm, x)
    # Dividing by the endpoint makes `quantile(d, 1) == x[end]` true by construction.
    return InterpolatedDistribution(x, y, norm, cdf ./ last(cdf))
end

"""
    normalizer(d) -> Real

Trapezoid integral of the tabulated (unnormalized) density.
"""
normalizer(d::InterpolatedDistribution) = d.norm

Base.minimum(d::InterpolatedDistribution) = first(d.x)
Base.maximum(d::InterpolatedDistribution) = last(d.x)
Base.eltype(d::InterpolatedDistribution) = promote_type(eltype(d.y), typeof(d.norm))

Distributions.insupport(d::InterpolatedDistribution, v::Real) = minimum(d) <= v <= maximum(d)

function Distributions.logpdf(d::InterpolatedDistribution, v::Real)
    return log(_interp(v, d.x, d.y; left = 0, right = 0) / d.norm)
end

Distributions.cdf(d::InterpolatedDistribution, v::Real) = _interp(v, d.x, d.cdf_grid)

Distributions.quantile(d::InterpolatedDistribution, q::Real) = _interp(q, d.cdf_grid, d.x)

function Random.rand(rng::AbstractRNG, d::InterpolatedDistribution)
    return quantile(d, rand(rng))
end
