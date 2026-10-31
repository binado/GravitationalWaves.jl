# GWMatchedFiltering

Whitened frequency-domain inner products, SNRs and fixed-PSD Gaussian likelihoods on plain
arrays. Data are `(frequency, channel)` complex arrays in one-sided physical Fourier
normalization; the noise enters through a whitening factor `W` and the frequency spacing `df`.

```julia
using GWMatchedFiltering, GWNoise
f = 0:0.25:1024
S = band_psd(TabulatedPSD(psd_path(:AplusDesign), ' '), f; fmin = 20, fmax = 512)
W = whitening_factor(reshape(S, :, 1))
h = ones(ComplexF64, length(f), 1)
data = copy(h)
log_likelihood_ratio(h, data, W, step(f))
```

## Whitening factor

`W = L⁻¹` with `L L† = S` per bin, zero on inactive bins, in one of two shapes:

- `(nf, nc)` real, `1/√S`, for independent channels (a vector for one channel);
- `(nc, nc, nf)` lower triangular, for correlated channels.

`GWNoise.whitening_factor` builds it from a PSD, validating the input once; no inverse
covariance is ever formed. The inner product is `⟨a, b⟩ = 4Δf Σₖ (W aₖ)† (W bₖ)`, conjugating
the first argument. Statistics read the shape of `W` through `filter_structure`, which returns
`DiagonalFilter()` or `BlockFilter()`; `whiten(d, W)` and `whiten!(out, d, W)` are the only
structure-dependent operations.

## API

`complex_inner_product(a, b, W, df)`, `noise_weighted_inner_product`, `optimal_snr_squared(h, W, df)`,
`optimal_snr`, `matched_filter_snr(template, data, W, df)`, `residual_squared(data, template, W, df)`,
`likelihood_statistics` and `log_likelihood_ratio`. Kernels check no shapes or values; `a` and `b`
must have the same shape. Matched-filter SNR assumes a finite nonzero template, so run
`check_template(template, W, df)` explicitly on the host.

## Differentiation and performance

ForwardDiff and Mooncake acceptance lives in the package’s `test/test_acceptance_*.jl` files; only strain
parameters are differentiated, `W` and `df` are fixed. Enzyme and Reactant are not supported
(see `AGENTS.md` for the pitfalls found). The kernels are written as array expressions so they
trace and differentiate, which allocates: for three channels they are roughly 4x (diagonal) to
6-16x (blocks) slower than hand-written host loops. A host fast path is a recorded TODO in
`src/whitened.jl`. Layout measurements live in `benchmarks/GWMatchedFiltering`.

## Local acceptance tests

The acceptance items check joint metrics, waveform/projection/likelihood composition, metadata and reconstruction boundaries, single-channel vectors and fixed Planck windows.

From the repository root, run `just test-integration-inference`. These `:integration` items are excluded
from default `Pkg.test()` and CI. See the [selection and SNR conventions](../../README.md#development-and-tests).
