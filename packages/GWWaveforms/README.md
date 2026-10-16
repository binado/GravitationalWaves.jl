# GWWaveforms

Pure-Julia frequency-domain compact-binary waveforms with the interface of
[`LALWaveforms`](../LALWaveforms): `WaveformSettings`, `hphc_fd`, and [bilby](https://bilby-dev.github.io/bilby/)'s
parameter names, units and mass conversions. The code is generic over `Real`, so it differentiates with
[ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl), [Enzyme](https://github.com/EnzymeAD/Enzyme.jl) and
[Mooncake](https://github.com/chalk-lab/Mooncake.jl), needs no LAL install, and the compiler can specialise it.

| approximant  | model |
|--------------|-------|
| `TaylorF2`   | post-Newtonian stationary-phase inspiral, aligned spins, optional SPA amplitude corrections |
| `IMRPhenomD` | aligned-spin inspiral-merger-ringdown (Husa et al. 2016, Khan et al. 2016) |

```julia
using GWWaveforms

settings = WaveformSettings("IMRPhenomD"; reference_frequency = 50.0, minimum_frequency = 20.0)
f = range(0, 1024; step = 0.25)
params = (
    chirp_mass = 28.0, mass_ratio = 0.8,   # or mass_1 and mass_2, detector-frame solar masses
    luminosity_distance = 500.0,           # Mpc
    theta_jn = 0.4, phase = 1.3,
    a_1 = 0.5, a_2 = 0.3, tilt_1 = 0.0, tilt_2 = Float64(π),   # (anti-)aligned spins only
)
(; plus, cross) = hphc_fd(f, params, settings)
```

The element type of the output follows the parameters, so dual numbers pass through:

```julia
using ForwardDiff

loss(x) = sum(abs2, hphc_fd(f, merge(params, (; chirp_mass = x[1], mass_ratio = x[2])), settings).plus)
ForwardDiff.gradient(loss, [28.0, 0.8])
```

## Scope

- Frequency domain, uniform `AbstractRange` grids starting at 0, detector-frame masses, no tidal parameters.
- Both models are aligned-spin and quasi-circular: spins must be (anti-)aligned (`tilt` 0 or π, or zero spin) and
  `eccentricity` zero, otherwise `ArgumentError` is thrown (LAL rejects such input too). Neither model has higher
  modes, so `mode_array` is rejected.
- `pn_spin_order`, `pn_phase_order` and `pn_amplitude_order` act as in LAL's TaylorF2 (`-1` selects the highest
  order); IMRPhenomD ignores them apart from the start-frequency shift of `pn_amplitude_order`.

## Fidelity

The ports follow LALSimulation (lalsuite 7.26) line by line, including the operation order, LAL's `UsefulPowers`
and the natural cubic spline that GSL builds over the QNM tables. Against `LALWaveforms` the largest deviation
relative to the largest sample is about 1e-14 to 4e-13 for IMRPhenomD and 1e-13 to 7e-12 for TaylorF2 (limited by
`cbrt`/`log`/`sin` ulp differences on phases of order 1e5 rad). Where LAL and a paper or ripple disagree, LAL wins.

## Tests

The comparison with LAL needs the LAL libraries (see the [repository README](../../README.md)):

```sh
export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)
julia --project=packages/GWWaveforms -e 'using Pkg; Pkg.test()'
```

Without them the comparison warns and is skipped (with `CI=true` it fails). The AD tests differentiate a scalar loss
with all three backends through DifferentiationInterface and compare against finite differences and each other.

## License

GPL-2.0-or-later: the models are ports of LALSimulation.
