# GWSurrogates

Differentiable Julia ports of the NRSur7dq4 / NRSur7dq4v2 numerical-relativity surrogates from
[`gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate), designed to work with
ForwardDiff, Enzyme and Mooncake.

**Status: kernels only.** The package currently holds the building blocks; the surrogate model
itself (fits, dynamics, waveform assembly) follows. The model is time-domain, so the first
release returns time-domain modes on a fixed time grid.

| Kernel | Notes |
|---|---|
| `NaturalCubicSpline`, `spline_eval` | Natural cubic spline for real or complex data, matching `gwsurrogate`'s `_spline_interp.cpp`. |
| `wigner_D`, `wigner_D_element` | Closed-form Wigner-D for `l <= 5` in the spinor entries `a = w + iz`, `b = y + ix`. No square roots, no branch on the rotor. |
| `guarded_sqrt`, `rescale_to_norm` | Norm and rescaling that stay finite (and differentiable) at zero spin. |

`SphericalFunctions.jl` and `DataInterpolations.jl` are test oracles only.

## Differentiation notes

- Exactly-zero spin is not a regular point: the renormalisation target `|chi|` has a cone-shaped
  kink there. `rescale_to_norm` returns a zero vector unchanged instead of `NaN`.
- Model data will be stored in flat integer-indexed vectors and read with plain loops: Enzyme
  fails on `Dict` lookups with runtime keys and on closures over arrays loaded from the model.

## Data

The model files are the `NRSur7dq4.h5` and `NRSur7dq4v2.h5` data of the Zenodo record
[10.5281/zenodo.22257361](https://doi.org/10.5281/zenodo.22257361), licensed CC-BY-4.0. Cite the
record and Field et al., [JOSS 10(107), 7073](https://joss.theoj.org/papers/10.21105/joss.07073)
if you use them.

## License

MIT; see [LICENSE](LICENSE).
