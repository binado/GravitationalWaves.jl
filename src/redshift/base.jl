using BackgroundCosmology: cosmology, distance_and_volume_grid
using Distributions
using Random

"""
    RedshiftDistribution(; params, source_frame_distribution, minimum_redshift = 0.0,
                         maximum_redshift = 10.0, n_grid = 1000)

Detector-frame redshift distribution induced by a source-frame merger-rate density
`ψ(z)`,

    p(z | θ) ∝ ψ(z) / (1 + z) · dV_c/dz,

tabulated on a linear grid and wrapped in an [`InterpolatedDistribution`](@ref), which owns
the normalization, `logpdf`, `cdf`, `quantile` and `rand`.

`params` is a `NamedTuple` holding the cosmology keys of BackgroundCosmology (`H0`, `Ωm`,
optionally `w0`, `wa`) plus whatever `source_frame_distribution(z, params)` reads; the
modified-propagation keys `xi_0` and `xi_n` are used by [`distance_ratio`](@ref) when
present. Note that `minimum_redshift = 0` puts a zero-density node at the left edge, since
`dV_c/dz` vanishes at `z = 0`.
"""
struct RedshiftDistribution{
        R, P <: NamedTuple, D <: InterpolatedDistribution, V <: AbstractVector,
    } <: ContinuousUnivariateDistribution
    merger_rate_fn::R
    params::P
    dist::D
    luminosity_distance_grid::V
    differential_comoving_volume_grid::V
end

function RedshiftDistribution(;
        params::NamedTuple,
        source_frame_distribution,
        minimum_redshift::Real = 0.0,
        maximum_redshift::Real = 10.0,
        n_grid::Integer = 1000
    )
    z = collect(range(minimum_redshift, maximum_redshift, n_grid))
    grids = distance_and_volume_grid(cosmology(params), z)
    dVc = grids.differential_comoving_volume
    rate = [source_frame_distribution(zi, params) for zi in z]
    y = @. rate / (1 + z) * dVc
    return RedshiftDistribution(
        source_frame_distribution, params, InterpolatedDistribution(z, y),
        grids.luminosity_distance, dVc
    )
end

normalizer(d::RedshiftDistribution) = normalizer(d.dist)
Base.minimum(d::RedshiftDistribution) = minimum(d.dist)
Base.maximum(d::RedshiftDistribution) = maximum(d.dist)
Base.eltype(d::RedshiftDistribution) = eltype(d.dist)
Distributions.insupport(d::RedshiftDistribution, v::Real) = insupport(d.dist, v)
Distributions.logpdf(d::RedshiftDistribution, v::Real) = logpdf(d.dist, v)
Distributions.cdf(d::RedshiftDistribution, v::Real) = cdf(d.dist, v)
Distributions.quantile(d::RedshiftDistribution, q::Real) = quantile(d.dist, q)
Random.rand(rng::AbstractRNG, d::RedshiftDistribution) = rand(rng, d.dist)

"""
    merger_rate(d, z)

Source-frame merger rate `ψ(z)` evaluated pointwise at the hyperparameters `d` was built
with (not interpolated from the grid). For a time-delayed distribution this is the delayed
rate.
"""
merger_rate(d::RedshiftDistribution, z::Real) = d.merger_rate_fn(z, d.params)

# The user-supplied shape: the rate itself, except for a time-delayed rate, where it is
# the formation rate.
_source_frame(f) = f

"""
    source_frame_distribution(d, z, params)

The supplied source-frame rate shape `ψ(z)` evaluated with `params` (the formation rate
for a time-delayed distribution).
"""
source_frame_distribution(d::RedshiftDistribution, z::Real, params) = _source_frame(d.merger_rate_fn)(z, params)

"""Redshift grid the table was built on (the interpolant abscissa)."""
redshift_grid(d::RedshiftDistribution) = d.dist.x

"""Number of grid nodes."""
n_grid(d::RedshiftDistribution) = length(d.dist.x)

"""
    total_merger_rate(d) -> Real

Total merger rate in mergers per second, `1e-9 ∫ ψ(z)/(1+z) dV_c/dz dz / yr`, with `ψ` in
Gpc⁻³ yr⁻¹.
"""
total_merger_rate(d::RedshiftDistribution) = 1.0e-9 * normalizer(d) / JULIAN_YEAR_SEC

"""Luminosity distance in Mpc at `z`, clamped outside the grid."""
luminosity_distance(d::RedshiftDistribution, z::Real) = _interp(z, d.dist.x, d.luminosity_distance_grid)

"""
    distance_ratio(d, z)

GW-to-EM luminosity-distance ratio `Ξ(z)`: from `xi_0` and `xi_n` when `params` has them,
else 1.
"""
function distance_ratio(d::RedshiftDistribution, z::Real)
    if haskey(d.params, :xi_0)
        return exp(_log_gw_em_ratio(z, d.params.xi_0, d.params.xi_n))
    end
    return one(float(z))
end

"""Gravitational-wave luminosity distance in Mpc: [`luminosity_distance`](@ref) times [`distance_ratio`](@ref)."""
gw_luminosity_distance(d::RedshiftDistribution, z::Real) = luminosity_distance(d, z) * distance_ratio(d, z)

"""All-sky `dV_c/dz` in Mpc³ at `z`, clamped outside the grid."""
differential_comoving_volume(d::RedshiftDistribution, z::Real) = _interp(z, d.dist.x, d.differential_comoving_volume_grid)
