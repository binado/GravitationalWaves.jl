# GWSurrogates

Differentiable Julia ports of the NRSur7dq4 / NRSur7dq4v2 numerical-relativity surrogates from
[`gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate), designed to work with
ForwardDiff, Enzyme and Mooncake.

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
out.t, out.quat, out.orbphase, out.chiA, out.chiB   # coprecessing frame, sampled at output_times(dyn)
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
| `t_ref_from_omega` | Reference time from an orbital frequency by linear interpolation between nodes. Only C⁰: its derivative jumps at every node crossing, so fix `t_ref` for gradient-based work. |
| `CoorbitalSurrogate`, `load_coorbital`, `coorbital_modes` | Coorbital-frame modes from node fits and the empirical-interpolation basis; v2 blends two time subdomains with stored masks. |

`SphericalFunctions.jl` and `DataInterpolations.jl` are test oracles only.

## Differentiation notes

- Exactly-zero spin is not a regular point: the renormalisation target `|chi|` has a cone-shaped
  kink there. `rescale_to_norm` returns a zero vector unchanged instead of `NaN`.
- Model data is stored in flat integer-indexed vectors (CSR-style fit sets) and read with plain
  loops: Enzyme fails on `Dict` lookups with runtime keys and on closures over arrays loaded from
  the model. Enzyme and Mooncake acceptance of the coorbital surrogate is still to be checked.

## Differences from gwsurrogate

`gwsurrogate`'s C++ spline (`_spline_interp.cpp`, `hunt`) picks the wrong interval for an isolated
evaluation point in the last interval of a short knot window: it evaluates the cubic of the
interval before it. In the dynamics surrogate this changes the interpolated ODE derivative only for
times in the last interval of the grid (`t` in the final ~2 M, after merger). This package uses the
correct natural spline there, so `dynamics` with a `t_ref` in that range differs from
`gwsurrogate`; everywhere else it agrees to about 1e-14.

## Data

The model files are the `NRSur7dq4.h5` and `NRSur7dq4v2.h5` data of the Zenodo record
[10.5281/zenodo.22257361](https://doi.org/10.5281/zenodo.22257361), licensed CC-BY-4.0. Cite the
record and Field et al., [JOSS 10(107), 7073](https://joss.theoj.org/papers/10.21105/joss.07073)
if you use them.

## License

MIT; see [LICENSE](LICENSE).
