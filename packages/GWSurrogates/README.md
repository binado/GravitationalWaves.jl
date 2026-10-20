# GWSurrogates

Differentiable Julia ports of the NRSur7dq4 / NRSur7dq4v2 numerical-relativity surrogates from
[`gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate), designed to work with
ForwardDiff and Mooncake.

**Status: coorbital waveform and dynamics surrogates.** The package holds the kernels, the
coorbital-frame waveform surrogate and the dynamics surrogate of NRSur7dq4 and NRSur7dq4v2; the
inertial-frame assembly and the public waveform generator follow. The model is time-domain, so the first release returns time-domain modes on a
fixed time grid.

```julia
using GWSurrogates
sur = GWSurrogates.load_coorbital(:NRSur7dq4)   # downloads and caches the data on first use
# chiA, chiB: coorbital-frame spins, (length(sur.t), 3)
h = GWSurrogates.coorbital_modes(sur, q, chiA, chiB)   # (n_times, n_modes), (2,-2), …, (4,4)
```

```julia
dyn = GWSurrogates.load_dynamics(:NRSur7dq4)
# Evolve from the reference time t_ref (default: first node) to the full time grid.
out = GWSurrogates.dynamics(dyn, q, chiA0, chiB0; t_ref = -100.0)
out.quat, out.orbphase, out.chiA, out.chiB   # coprecessing frame, sampled at output_times(dyn)
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
| `DynamicsSurrogate`, `load_dynamics`, `dynamics` | The 11-variable precession ODE: RK4 start-up then variable-step AB4 on the fixed node grid, forward and backward from `t_ref`, with per-step renormalization. Allocation-light (`SVector` state) and about 0.6 ms per trajectory. |
| `PrecessingSurrogate`, `load_surrogate`, `inertial_modes`, `mode_sum` | The full pipeline: dynamics interpolated onto the coorbital grid, coorbital modes, Wigner-D rotation into the inertial frame, spline resampling and the `₋₂Y_ℓm` sum. About 15 ms for 2000 time samples. |
| `PrecessingGenerator`, `hphc_td` | SI-unit `h₊`, `h×` from total mass, distance, mass ratio, spin vectors, inclination and reference phase. |
| `t_ref_from_omega` | Reference time from an orbital frequency by linear interpolation between nodes. Only C⁰: its derivative jumps at every node crossing, so fix `t_ref` for gradient-based work. |
| `CoorbitalSurrogate`, `load_coorbital`, `coorbital_modes` | Coorbital-frame modes from node fits and the empirical-interpolation basis; v2 blends two time subdomains with stored masks. |

`SphericalFunctions.jl` and `DataInterpolations.jl` are test oracles only.

## Differentiation

ForwardDiff and Mooncake both differentiate the kernels,
the coorbital surrogate, the dynamics integration and the full strain for both models, and agree
with finite differences (`just test-integration-surrogates`, a local suite in
`integration/GWSurrogates`; CI does not run it). Reactant is out of scope.

Cost of the gradient of a scalar functional of `hphc_td` (`ell_max = 4`, about 2600 samples at
2048 Hz, 8 parameters: mass, mass ratio, six spin/orientation parameters), Apple M-series CPU,
Julia 1.12:

| | primal | ForwardDiff | Mooncake |
|---|---|---|---|
| NRSur7dq4 | 7 ms | 70 ms | 133 ms |
| NRSur7dq4v2 | 9 ms | 76 ms | 214 ms |
| first-call compile | – | ≈ 0.1 s | ≈ 42 s |

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
