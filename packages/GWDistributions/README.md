# GWDistributions.jl

[![CI](https://github.com/binado/GravitationalWaves.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/binado/GravitationalWaves.jl/actions/workflows/CI.yml)

Population-distribution building blocks for gravitational-wave source modeling, ported from the
Python/NumPyro `astrogwb.distributions` and built on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl):

- `RedshiftDistribution`: detector-frame redshift density `ψ(z)/(1+z) · dV_c/dz` for any source-frame
  rate `ψ`, with `madau_dickinson_redshift_distribution` and
  `madau_dickinson_time_delayed_redshift_distribution` (mergers delayed from a formation rate by a
  delay-time distribution) as ready-made factories;
- `InterpolatedDistribution`: tabulated density with trapezoid normalizer, `cdf` and `quantile`;
- `MaxOfTwoNormals`, `MaxOfTwoUniforms`, `PowerLawDelay` (smooth through `α = -1`) and
  `UniformCosine`, each with closed-form `cdf`/`quantile`;
- the default BBH mass law (`DefaultBBHPrimaryMass` / `DefaultBBHMassPair` with `planck_taper`).

Cosmology comes from [BackgroundCosmology.jl](https://github.com/binado/BackgroundCosmology.jl).
`params` is a `NamedTuple` holding its keys (`H0`, `Ωm`, optionally `w0`, `wa`) next to the rate
parameters (`γ`, `κ`, `z_peak`, optional `local_merger_rate`) and the optional modified-propagation
keys `xi_0`, `xi_n`, so a single tuple feeds both `cosmology(params)` and the distribution.

## Installation

```julia
using Pkg
Pkg.add([
    PackageSpec(url = "https://github.com/binado/BackgroundCosmology.jl"),
    PackageSpec(url = "https://github.com/binado/GravitationalWaves.jl", subdir = "packages/GWDistributions"),
])
```

## Example

```julia
using Distributions
using GWDistributions

params = (H0 = 67.66, Ωm = 0.30966, γ = 2.7, κ = 3.0, z_peak = 2.0, local_merger_rate = 17.0)  # R₀ in Gpc⁻³ yr⁻¹

d = madau_dickinson_redshift_distribution(; params)
logpdf(d, 0.5)           # detector-frame density at z = 0.5
rand(d, 10)              # 10 redshift samples
total_merger_rate(d)     # mergers per second
luminosity_distance(d, 0.5)

# Mergers delayed from the same formation rate by p(τ) ∝ 1/τ on [0.02, 13] Gyr
td = madau_dickinson_time_delayed_redshift_distribution(;
    params, time_delay_distribution = PowerLawDelay(-1.0, 0.02, 13.0),
)
merger_rate(td, 0.0)     # == local_merger_rate
```

Everything is differentiable with ForwardDiff: quadratures have a fixed order and the cosmology
inversion uses a fixed iteration count.

## Development

Run these commands from `packages/GWDistributions` in the monorepo:

```
just install-runic   # once
just fmt             # format with Runic
just test            # Pkg.test() on Julia 1.12
```

## License

MIT; see [LICENSE](LICENSE).
