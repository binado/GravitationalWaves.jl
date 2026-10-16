# GWDistributions.jl

[![CI](https://github.com/binado/GWDistributions.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/binado/GWDistributions.jl/actions/workflows/CI.yml)

Population-distribution building blocks for gravitational-wave source modeling: interpolated 1-D
distributions with cached normalizers (`Interpolated1DDistribution`), the DEFAULT BBH mass law
(`DefaultBBHPrimaryMass` / `DefaultBBHMassPair` with `planck_taper`), and the Madau–Dickinson
source-frame merger rate (`MadauDickinsonSourceFrame`) composed into a detector-frame redshift
distribution (`RedshiftInterpolatedDistribution`) whose `normalizer` is the total merger rate in
events/sec once the local rate `R₀` is included.

The package is cosmology-independent by design: redshift distributions take solid-angle-integrated
`dV_c/dz` **arrays**, not cosmology types. Pair it with
[BackgroundCosmology.jl](https://github.com/binado/BackgroundCosmology.jl) (not a dependency) for
real cosmology.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/binado/GWDistributions.jl")
```

## Example

A Madau–Dickinson redshift distribution from a volume grid, with `logpdf` and `rand`:

```julia
using GWDistributions
using BackgroundCosmology   # or any other source of dV_c/dz

c = cosmology(; H0 = 67.0, Ωm = 0.315)
z = collect(LinRange(0.0, 2.0, 201))
dVc = distance_and_volume_grid(c, z).differential_comoving_volume   # all-sky, Mpc³ per unit z

sf = MadauDickinsonSourceFrame(; γ = 2.7, κ = 3.0, zpeak = 2.5, R₀ = 161.0)  # R₀ in Gpc⁻³ yr⁻¹
d = RedshiftInterpolatedDistribution(sf, dVc, z)

logpdf(d, 0.5)    # detector-frame density at z = 0.5
rand(d, 10)       # 10 redshift samples
normalizer(d)     # total detector-frame merger rate, events/sec
```

## Development

```
just install-runic   # once
just fmt             # format with Runic
just test            # Pkg.test() on Julia 1.12
```
