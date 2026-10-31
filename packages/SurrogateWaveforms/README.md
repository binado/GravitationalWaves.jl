# SurrogateWaveforms

Differentiable Julia ports of the NRSur7dq4 / NRSur7dq4v2 numerical-relativity surrogates from
[`gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate), designed to work with
ForwardDiff and Mooncake.

**Status: coorbital waveform and dynamics surrogates.** The package holds the kernels, the
coorbital-frame waveform surrogate and the dynamics surrogate of NRSur7dq4 and NRSur7dq4v2; the
inertial-frame assembly and the public waveform generator follow. The model is time-domain, so the first release returns time-domain modes on a
fixed time grid.

```julia
using SurrogateWaveforms
sur = SurrogateWaveforms.load_coorbital(:NRSur7dq4)   # downloads and caches the data on first use
# chiA, chiB: coorbital-frame spins, (length(sur.t), 3)
h = SurrogateWaveforms.coorbital_modes(sur, q, chiA, chiB)   # (n_times, n_modes), (2,-2), …, (4,4)
```

```julia
dyn = SurrogateWaveforms.load_dynamics(:NRSur7dq4)
# Evolve from the reference time t_ref (default: first node) to the full time grid.
out = SurrogateWaveforms.dynamics(dyn, q, chiA0, chiB0; t_ref = -100.0)
out.t, out.quat, out.orbphase, out.chiA, out.chiB   # coprecessing frame, sampled at out.t = output_times(dyn)
```

Set `GWSURROGATES_DATA_DIR` to a directory holding `NRSur7dq4.h5` / `NRSur7dq4v2.h5` to skip the
download. Tests compare against `gwsurrogate`; generate the fixtures with
`uv run --project python/surrogates python python/generate_surrogates.py --data-dir DIR --download`.

| Kernel | Notes |
|---|---|
| `NaturalCubicSpline`, `spline_eval` | Natural cubic spline for real or complex data, matching `gwsurrogate`'s `_spline_interp.cpp`. |
| `wigner_D`, `wigner_D_element` | Closed-form Wigner-D for `l <= 5` in the spinor entries `a = w + iz`, `b = y + ix`. No square roots, no branch on the rotor. |
| `guarded_sqrt`, `rescale_to_norm` | Norm and rescaling that stay finite (and differentiable) at zero spin. |
| `FitSet`, `fit_powers!`, `eval_fit`, `fit_params` | Parametric fits (sum of `coef * prod x_j^n_j` over 7 parameters) in flat, integer-indexed storage. |
| `DynamicsSurrogate`, `load_dynamics`, `dynamics` | The 11-variable precession ODE: RK4 start-up then variable-step AB4 on the fixed node grid, forward and backward from `t_ref`, with per-step renormalization. Allocation-light (`SVector` state) and about 0.5 ms per trajectory. |
| `PrecessingSurrogate`, `load_surrogate`, `inertial_modes`, `mode_sum` | The full pipeline: dynamics interpolated onto the coorbital grid, coorbital modes, Wigner-D rotation into the inertial frame, spline resampling and the `₋₂Y_ℓm` sum. About 2 ms (NRSur7dq4) and 5 ms (NRSur7dq4v2, up to `ℓ = 5`) for 2000 time samples. |
| `PrecessingGenerator`, `hphc_td` | SI-unit `h₊`, `h×` from total mass, distance, mass ratio, spin vectors, inclination and reference phase. |
| `t_ref_from_omega` | Reference time from an orbital frequency by linear interpolation between nodes. Only C⁰: its derivative jumps at every node crossing, so fix `t_ref` for gradient-based work. |
| `CoorbitalSurrogate`, `load_coorbital`, `coorbital_modes` | Coorbital-frame modes from node fits and the empirical-interpolation basis; v2 blends two time subdomains with stored masks. |

`SphericalFunctions.jl` and `DataInterpolations.jl` are test oracles only.

## Differentiation

ForwardDiff and Mooncake both differentiate the kernels,
the coorbital surrogate, the dynamics integration and the full strain for both models, and agree
with finite differences (`just test-integration-surrogates`, a local suite in
the package’s `test/test_acceptance_*.jl` files; CI does not run it). Enzyme and Reactant are out of scope.

Cost of the gradient of a scalar functional of `hphc_td` (`ell_max = 4`, about 2600 samples at
2048 Hz, 8 parameters: mass, mass ratio, six spin/orientation parameters), Apple M-series CPU,
Julia 1.12:

| | primal | ForwardDiff | Mooncake |
|---|---|---|---|
| NRSur7dq4 | 2.2 ms | 22 ms | 83 ms |
| NRSur7dq4v2 | 3.4 ms | 25 ms | 130 ms |
| first-call compile | ≈ 1.5 s | ≈ 6 s | ≈ 57 s |

The compile row is for NRSur7dq4 in a fresh session; the Mooncake figure includes its one-time
warm-up, and NRSur7dq4v2 then needs about 8 s more. The Wigner-D rotation of each `ℓ` block is
generated code (`@generated`, unrolled sums with literal indices): it is 2.2 times faster than
loops with runtime bounds and bit-identical to them, but it adds roughly 20 to 30 s to the
Mooncake compile time (adding it lengthened `just test-integration-surrogates` by about 30 s).

With 8 parameters ForwardDiff is competitive; reverse mode pays off with more parameters (e.g. adding the initial frame or several `ell_max` studies).

Rules the implementation follows (each was a real failure while porting):

- Model data lives in flat integer-indexed vectors read with plain loops; no `Dict`, no
  closures or generators over arrays loaded from the model.
- Exactly-zero spin is not a regular point (the renormalisation target `|chi|` has a cone-shaped
  kink there); `rescale_to_norm` returns a zero vector unchanged instead of `NaN`.
- No branches on differentiated values inside a recurrence (the Thomas sweep of the spline starts
  from zeroed buffers instead of special-casing its first row).
- Derive the element type of a spin tuple with `element_type`: `eltype` of a tuple that mixes dual
  numbers and floats is an abstract union and slows everything down by orders of magnitude.
- `t_ref_from_omega` (used when you pass `f_ref`) is only C⁰, so pass `t_ref` when differentiating.

## Differences from gwsurrogate

`gwsurrogate`'s C++ spline (`_spline_interp.cpp`, `hunt`) picks the wrong interval for an isolated
evaluation point in the last interval of a short knot window: it evaluates the cubic of the
interval before it. In the dynamics surrogate this changes the interpolated ODE derivative only for
times in the last interval of the grid (`t` in the final ~2 M, after merger). This package uses the
correct natural spline there, so `dynamics` with a `t_ref` in that range differs from
`gwsurrogate`; everywhere else it agrees to about 1e-14. The same bug affects `gwsurrogate`'s
final resampling when an output sample falls in the last interval of the model grid (`t` within
about 0.2 M of the end), so the last output sample can differ from `gwsurrogate` by up to about 1e-6 of
the peak amplitude when the output grid ends exactly at the end of the model.

## Data

The model files are the `NRSur7dq4.h5` and `NRSur7dq4v2.h5` data of the Zenodo record
[10.5281/zenodo.22257361](https://doi.org/10.5281/zenodo.22257361), licensed CC-BY-4.0. Cite the
record and Field et al., [JOSS 10(107), 7073](https://joss.theoj.org/papers/10.21105/joss.07073)
if you use them.

## License

MIT; see [LICENSE](LICENSE).


### Shared waveform interface

`PrecessingGenerator <: AbstractWaveforms.HigherModesWaveformGenerator` implements
`waveform_modes_td(generator, params; dt=nothing, times=nothing)`, returning seconds, an
ordered tuple of signed `(l,m)` labels, and a `(sample, mode)` complex matrix in physical
strain. Columns run from `l=2` through `ell_max`, with `m=-l:l`. These are the inertial
modes in the LAL source frame at the reference epoch; `init_orbphase` and `init_quat`
retain their existing frame interpretation. Modes exclude both `inclination` and `phi_ref`.
The observer convention is `(theta, phi) = (inclination, π/2 - phi_ref)`.

`hphc_td` uses shared harmonic assembly and requires `dt` or a uniform `times` range and returns a two-channel TimeSeries instead of `(t, hplus, hcross)`
result. `inertial_modes`, `mode_sum(modes, ell_max, inclination, phi_ref)`, and `sylm_m2`
remain available. The generator reports precessing spin support and differentiable evaluation;
its frequency accessors return `f_ref`, `f_low`, and `nothing` for the upper bound.
It does not provide FD modes or arbitrary mode-array selection.

Run `julia --project=packages/SurrogateWaveforms benchmarks/SurrogateWaveforms/assembly.jl` to measure
harmonic assembly separately from generation and compare the previous Wigner-D sum.

For assembled strain use `times(h)`, `h[:, :plus]`, and `h[:, :cross]` from GWSeries.
Lower-level `waveform_modes_td` and `inertial_modes` retain native irregular grids.

Coordinate-first `hphc_td(times, params, generator)` returns an ordinary matrix with
plus/cross columns and supports the same time coordinates as `waveform_modes_td`.
`hphc_td(prototype::TimeSeries, params, generator)` returns a TimeSeries on the
prototype's peak-relative time grid with channels `(:plus, :cross)`. The existing
`hphc_td(generator, params; dt=...)` and `times=uniform_range` conveniences still
return TimeSeries. Detector GPS arrival-time placement is separate from these calls.

## Local acceptance tests

Acceptance checks native kernels, both models’ coorbital modes, dynamics, fixed-reference-time waveform gradients and SI strain, plus windowed and padded TD SNR. It requires model data from Scratch or `GWSURROGATES_DATA_DIR`, and needs no Python. Shared harmonic assembly AD belongs to AbstractWaveforms. The live upstream comparison belongs to [SurrogatePyCallWaveforms](../SurrogatePyCallWaveforms/README.md).

From the repository root, run `just test-integration-surrogates`. These `:integration` items are excluded
from default `Pkg.test()` and CI. See the [selection and SNR conventions](../../README.md#development-and-tests).
