# FrequencyDomainWaveforms

Pure-Julia frequency-domain compact-binary waveforms with the interface of
[`LALWaveforms`](../LALWaveforms): the generic `hphc_fd(f, params, generator)` of [`AbstractWaveforms`](../AbstractWaveforms), and [bilby](https://bilby-dev.github.io/bilby/)'s
parameter names, units and mass conversions. The code is generic over `Real`, so it differentiates with
[ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and
[Mooncake](https://github.com/chalk-lab/Mooncake.jl), needs no LAL install, and the compiler can specialise it.

| generator    | model |
|--------------|-------|
| `TaylorF2()` | post-Newtonian stationary-phase inspiral, aligned spins, optional SPA amplitude corrections |
| `IMRPhenomD()` | aligned-spin inspiral-merger-ringdown (Husa et al. 2016, Khan et al. 2016) |
| `IMRPhenomXAS()` | default aligned-spin quadrupole inspiral-merger-ringdown (Pratten et al. 2020) |

```julia
using FrequencyDomainWaveforms

generator = IMRPhenomD(; reference_frequency = 50.0, minimum_frequency = 20.0)
f = range(0, 1024; step = 0.25)
params = (
    chirp_mass = 28.0, mass_ratio = 0.8,   # or mass_1 and mass_2, detector-frame solar masses
    luminosity_distance = 500.0,           # Mpc
    theta_jn = 0.4, phase = 1.3,
    a_1 = 0.5, a_2 = 0.3, tilt_1 = 0.0, tilt_2 = Float64(π),   # (anti-)aligned spins only
)
h = hphc_fd(f, params, generator)  # ordinary matrix: plus/cross columns
(; amplitude, phase) = amplitude_phase_fd(f, params, generator)
```

Each model is its own generator type with only the options it uses; the generators and `hphc_fd` are the
ones `LALWaveforms` shares, so code can switch between the native and the LAL-backed implementation.

All three generators subtype `AlignedSpinQuadrupolarWaveformGenerator`. Their
`amplitude_phase_fd` methods return nonnegative amplitude (including distance normalization)
and analytic, unwrapped phase in radians. The inclination-independent strain is
`amplitude .* exp.(im .* phase)`; `hphc_fd` applies the shared inclination factors.
Both vectors are zero at DC and outside waveform support, with phase set to zero wherever
amplitude vanishes. TaylorF2's signed LAL amplitude is absorbed into the phase by a π shift.
Changing the inclination does not change amplitude or phase.

Each model implements `amplitude_phase_fd` directly. The Phenom models prepare their coefficients
once, then a plain loop evaluates one scalar kernel per frequency and fills the amplitude and
phase vectors, which are returned in a named tuple. The method trusts its inputs; run
`validate_parameters` on the host first.

The element type of the output follows the parameters, so dual numbers pass through:

```julia
using ForwardDiff

loss(x) = sum(abs2, hphc_fd(f, merge(params, (; chirp_mass = x[1], mass_ratio = x[2])), generator)[:, 1])
ForwardDiff.gradient(loss, [28.0, 0.8])
```

Parameter values are not checked by `hphc_fd` or `amplitude_phase_fd`. Call
`validate_parameters(f, params, generator)` once on the host: it throws `ArgumentError` for unknown
keys, spins that are not (anti-)aligned, a nonzero `eccentricity`, non-positive masses or distance,
and frequency ranges the model cannot produce, and returns `nothing` otherwise.

## Automatic differentiation

AD support is restricted to ForwardDiff and Mooncake. The frequency grid and generator
are fixed configuration; masses, spins, distance, inclination and phase may be differentiated.
Kernels use ordinary Julia loops, branches and small matrix solves.
IMRPhenomD uses component masses directly for mass differences and remnant fits,
keeping component-mass derivatives finite at exact equal masses.
`maximum_frequency` must be positive; intrinsic mass-dependent cutoffs are applied within
the configured frequency interval.

## Scope

- Frequency domain, uniform `AbstractRange` grids, detector-frame masses, no tidal parameters.
  All models evaluate the supplied frequencies directly and accept grids starting above zero,
  without a DC bin.
  Grid validation runs in `validate_parameters`, outside waveform evaluation.
- All models are aligned-spin and quasi-circular: spins must be (anti-)aligned (`tilt` 0 or π, or zero spin) and
  `eccentricity` zero, otherwise `validate_parameters` throws `ArgumentError` (LAL rejects such input too). None has higher
  modes, so there is no `mode_array` option.
- `TaylorF2` takes `pn_spin_order`, `pn_phase_order` and `pn_amplitude_order` as in LAL (`-1` selects the highest
  order); `IMRPhenomD` and `IMRPhenomXAS` only take `pn_amplitude_order`, which shifts the start frequency.
- `IMRPhenomXAS` implements LAL's default calibration versions. Its intrinsic cutoff is
  `0.3 / M_sec`, inclusive. Alternative calibration flags are not exposed.

### IMRPhenomXAS differentiation

The collocation systems use StaticArrays matrices and pivoted solves, differentiating both
the matrix and the right-hand side. The intermediate amplitude uses
a normalized polynomial coordinate to avoid large powers of nearby frequency nodes.

Mass differences are calculated directly from component masses, avoiding the singular derivative
of `sqrt(1 - 4eta)` at equal mass. The ISCO expression has a regular evaluation near zero final
spin. These choices do not remove coordinate singularities in conversions from other mass
parameters. Hard frequency cutoffs remain discontinuous: gradients describe the currently
selected support and do not differentiate a bin entering or leaving it.

## Fidelity

The TaylorF2 and IMRPhenomD ports follow LALSimulation (lalsuite 7.26) line by line, retaining the numerical kernels, LAL's `UsefulPowers`
and the natural cubic spline that GSL builds over the QNM tables. The public amplitude/phase
convention absorbs signed prefactors into phase; reconstruction can differ from LAL by
phase-rounding errors. Against `LALWaveforms` the largest deviation
relative to the largest sample is about 1e-14 to 4e-13 for IMRPhenomD and 1e-13 to 7e-12 for TaylorF2 (limited by
`cbrt`/`log`/`sin` ulp differences on phases of order 1e5 rad). Where LAL and a paper or ripple disagree, LAL wins.

IMRPhenomXAS fit tables and arithmetic were transcribed from
[Ripple revision 58c5302](https://github.com/tedwards2412/ripple/tree/58c53024718697e0430b666ef5dd40ff3c2e029d/src/ripplegw/waveforms),
with full complex polarizations and exact support compared against LALSimulation 7.26.
The comparison covers equal and unequal masses, ordered component labels, high spins,
high mass ratio, moving cutoffs and alternative frequency settings. LAL is the value oracle.
Across these cases, the largest complex-polarization error relative to the largest sample
is below `8e-12` (below `2e-13` for the tested ordinary stellar-mass binaries).

## Tests

The comparison with LAL needs the LAL libraries (see the [repository README](../../README.md)):

```sh
export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)
julia --project=packages/FrequencyDomainWaveforms -e 'using Pkg; Pkg.test()'
```

Without them the comparison warns and is skipped (with `CI=true` it fails). The default
AD tests differentiate a scalar loss with ForwardDiff through DifferentiationInterface
and compare against finite differences.

Run the separate integration suite from the repository root to check ForwardDiff
and Mooncake against finite differences and each other:

```sh
just test-integration-waveforms
```

This installs its own dependencies on first use and runs locally, outside CI. It
needs no LAL libraries.

## License

GPL-2.0-or-later: the models are ports of LALSimulation.
