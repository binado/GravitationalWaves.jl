using BackgroundCosmology: cosmology, lookback_time
using Distributions
using FastGaussQuadrature: gausslegendre

"""
    TimeDelayedMergerRate(formation_rate, delay, gl_rule, maximum_formation_redshift)

Callable `(z, params) -> R_m(z)`: the merger rate of binaries formed at rate
`formation_rate(z, params)` and merged after a delay `τ ~ delay` (Gyr). Being a convolution
in lookback time,

    R_m(z) ∝ ∫ p(τ) ψ(z_f(t_L(z) + τ)) dτ,

it is evaluated in the delay's probability variable `u = F(τ)` as
`∫₀^{F(τ_avail)} ψ(z_f(t_L(z) + F⁻¹(u))) du` with `τ_avail = t_L(z_cut) - t_L(z)`, using the
Gauss–Legendre rule `gl_rule` (from `gausslegendre`). The upper limit moves continuously
with the parameters (no mask), and the result is rescaled so `R_m(0) = ψ(0)`.
"""
struct TimeDelayedMergerRate{F, D <: UnivariateDistribution, W, T <: Real}
    formation_rate::F
    delay::D
    gl_rule::W
    maximum_formation_redshift::T
end

_source_frame(r::TimeDelayedMergerRate) = r.formation_rate

function _delayed_integral(r::TimeDelayedMergerRate, c, t_cut, z, params)
    t_merger = lookback_time(c, z)
    # Clip into the support so `cdf` saturates at 0 and 1.
    available = clamp(t_cut - t_merger, minimum(r.delay), maximum(r.delay))
    u, w = _gauss_legendre_nodes_weights(zero(available), cdf(r.delay, available), r.gl_rule)
    acc = zero(first(w) * r.formation_rate(zero(t_merger), params))
    for k in eachindex(u, w)
        z_f = _redshift_at_lookback_time(c, t_merger + quantile(r.delay, u[k]))
        acc += w[k] * r.formation_rate(z_f, params)
    end
    return acc
end

function (r::TimeDelayedMergerRate)(z::Real, params::NamedTuple)
    c = cosmology(params)
    t_cut = lookback_time(c, r.maximum_formation_redshift)
    local_rate = r.formation_rate(zero(float(z)), params)
    return _delayed_integral(r, c, t_cut, z, params) * local_rate /
        _delayed_integral(r, c, t_cut, zero(float(z)), params)
end

"""
    time_delayed_merger_rate(z, params; formation_rate, delay, n_delay_nodes = 48,
                             maximum_formation_redshift = 20.0)

Pointwise merger rate `R_m(z)` of `formation_rate(z, params)` delayed by `delay`; see
[`TimeDelayedMergerRate`](@ref). `params` must carry the cosmology keys.
"""
function time_delayed_merger_rate(
        z::Real, params::NamedTuple;
        formation_rate, delay, n_delay_nodes::Integer = 48,
        maximum_formation_redshift::Real = 20.0
    )
    rate = TimeDelayedMergerRate(
        formation_rate, delay, gausslegendre(n_delay_nodes), maximum_formation_redshift
    )
    return rate(z, params)
end

"""
    TimeDelayedRedshiftDistribution(; params, source_frame_distribution,
        time_delay_distribution, n_delay_nodes = 48, maximum_formation_redshift = 20.0,
        minimum_redshift = 0.0, maximum_redshift = 10.0, n_grid = 1000)

[`RedshiftDistribution`](@ref) whose merger rate is `source_frame_distribution` (read as
the *formation* rate `ψ`) delayed by `time_delay_distribution` (in Gyr, e.g.
`PowerLawDelay(-1, 0.02, 13)`); see [`TimeDelayedMergerRate`](@ref). The delayed rate is
rescaled to `R_m(0) = ψ(0)`, so `local_merger_rate` keeps meaning the merger rate today.
Formation is cut off above `maximum_formation_redshift`, independently of
`maximum_redshift`, which bounds only the merger redshift.
"""
const TimeDelayedRedshiftDistribution = RedshiftDistribution{<:TimeDelayedMergerRate}

function (::Type{T})(;
        params::NamedTuple,
        source_frame_distribution,
        time_delay_distribution,
        n_delay_nodes::Integer = 48,
        maximum_formation_redshift::Real = 20.0,
        kwargs...
    ) where {T <: RedshiftDistribution{<:TimeDelayedMergerRate}}
    rate = TimeDelayedMergerRate(
        source_frame_distribution, time_delay_distribution,
        gausslegendre(n_delay_nodes), maximum_formation_redshift
    )
    return RedshiftDistribution(; params, source_frame_distribution = rate, kwargs...)
end
