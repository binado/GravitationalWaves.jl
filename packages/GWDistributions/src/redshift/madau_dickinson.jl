using BackgroundCosmology: AbstractCosmology, comoving_volume_element

"""
    madau_dickinson_rate(z, γ, κ, z_peak, local_merger_rate = 1)

Madau–Dickinson source-frame rate density,

    ψ(z) = R₀ C (1+z)^γ / (1 + ((1+z)/(1+z_p))^(γ+κ)),   C = 1 + (1+z_p)^(-(γ+κ)),

in Gpc⁻³ yr⁻¹ for `R₀ = local_merger_rate`. With the default `R₀ = 1` it is the
dimensionless shape with `ψ(0) = 1`.
"""
function madau_dickinson_rate(z::Real, γ::Real, κ::Real, z_peak::Real, local_merger_rate::Real = 1.0)
    one_plus_z = 1 + z
    exponent = γ + κ
    normalization = 1 + (1 + z_peak)^(-exponent)
    return local_merger_rate *
        (normalization * one_plus_z^γ / (1 + (one_plus_z / (1 + z_peak))^exponent))
end

"""
    total_merger_rate(z, weights, merger_rate, c::AbstractCosmology) -> Real

Total merger rate in mergers per second by quadrature over the nodes `z` with `weights`,
`1e-9 Σⱼ wⱼ ψ(zⱼ)/(1+zⱼ) dV_c/dz(zⱼ) / yr`, where `merger_rate` is `ψ` at `z` in
Gpc⁻³ yr⁻¹. The same quantity as `total_merger_rate(::RedshiftDistribution)` with no
distribution or table.
"""
function total_merger_rate(
        z::AbstractVector{<:Real}, weights::AbstractVector{<:Real},
        merger_rate::AbstractVector{<:Real}, c::AbstractCosmology
    )
    volume = @. 4π * comoving_volume_element(c, z)
    return 1.0e-9 * sum(@. weights * merger_rate / (1 + z) * volume) / JULIAN_YEAR_SEC
end

# Adapts the `params` NamedTuple to `madau_dickinson_rate`.
function _madau_dickinson_source_frame(z, params)
    return madau_dickinson_rate(
        z, params.γ, params.κ, params.z_peak, get(params, :local_merger_rate, 1.0)
    )
end

"""
    madau_dickinson_redshift_distribution(; params, minimum_redshift = 0.0,
        maximum_redshift = 10.0, n_grid = 1000)

[`RedshiftDistribution`](@ref) for a Madau–Dickinson merger rate. `params` carries `γ`,
`κ`, `z_peak` and optionally `local_merger_rate`, on top of the cosmology keys.
"""
function madau_dickinson_redshift_distribution(; params::NamedTuple, kwargs...)
    return RedshiftDistribution(;
        params, source_frame_distribution = _madau_dickinson_source_frame, kwargs...
    )
end

const MadauDickinsonRedshiftDistribution = madau_dickinson_redshift_distribution

"""
    madau_dickinson_time_delayed_redshift_distribution(; params, time_delay_distribution,
        n_delay_nodes = 48, maximum_formation_redshift = 20.0, minimum_redshift = 0.0,
        maximum_redshift = 10.0, n_grid = 1000)

Mergers delayed from a Madau–Dickinson *formation* rate; `local_merger_rate` stays the
merger rate at `z = 0`. See [`TimeDelayedRedshiftDistribution`](@ref).
"""
function madau_dickinson_time_delayed_redshift_distribution(;
        params::NamedTuple, time_delay_distribution, kwargs...
    )
    return TimeDelayedRedshiftDistribution(;
        params, source_frame_distribution = _madau_dickinson_source_frame,
        time_delay_distribution, kwargs...
    )
end
