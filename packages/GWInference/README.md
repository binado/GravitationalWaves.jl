# GWInference

Joint frequency-domain noise metrics, inner products, and SNRs for Julia.
Runtime dependencies: `LinearAlgebra` and `StaticArrays`. Inputs are already projected detector
strains; no waveform generator, detector geometry, sampler, Python or LAL is
required.

```julia
using GWInference

f = 0.0:0.25:1024.0
channels = (:H1, :L1)
S = fill(1.0e-46, length(f), 2)  # one-sided PSDs, strain²/Hz
metric = FrequencyDomainMetric(f, channels, S; minimum_frequency=20.0, nsamples=8192)

h = FrequencyDomainStrain(f, channels, hcat(h1, l1))  # projected Fourier templates
d = FrequencyDomainStrain(f, channels, hcat(data_h1, data_l1))
optimal_snr(h, metric)
matched_filter_snr(h, d, metric)  # complex scalar at the supplied alignment
noise_weighted_inner_product(h, d, metric)  # real scalar, summed over channels
stats = likelihood_statistics(h, d, metric)
log_likelihood_ratio(h, d, metric)
```

`h1`, `l1`, `data_h1` and `data_l1` above are caller-supplied vectors. Raw matrices
can be passed directly to all statistics; a vector denotes a single channel.
Containers validate grid and channel identity at evaluation. Raw arrays implicitly
use the metric's ordering. Strain storage is retained without copying; metric
setup copies weights/factors. Treat prepared objects and their storage as fixed.

## Correlated channels

Supply spectral matrices shaped `(frequency, channel, channel)`:

```julia
S = zeros(ComplexF64, length(f), 3, 3)
for k in eachindex(f)
    S[k, :, :] = 1.0e-46 .* [1 0.2 0.2; 0.2 1 0.2; 0.2 0.2 1]
end
metric = FrequencyDomainMetric(f, (:ET1, :ET2, :ET3), S; minimum_frequency=20.0)
```

This is an illustrative correlated three-channel model, not an ET noise forecast.
Auto-spectra occupy the diagonal; cross-spectra must be supplied with Hermitian
symmetry. For interior bins define `S[k,i,j] = 2Δf E[n[k,i] conj(n[k,j])]`.
Stationarity makes the metric block-diagonal over frequencies, with detector
coupling inside each block. Independent noise is the diagonal specialization.

Setup caches inverse lower Cholesky factors `W[k] = L[k] \ I`, where
`S[k] = L[k] * L[k]'`. Contractions use `(W*a)' * (W*b)`; no global matrix or
inverse spectral covariance is formed. Active blocks must be positive definite.
Singular or indefinite blocks fail; there is no implicit jitter, pseudoinverse,
or deletion of correlations. A signal null stream alone does not imply singular
noise covariance.

## Numerical contract

- Physical Fourier strain is `rfft(strain)/sampling_frequency`.
- The complex contraction is `C(a,b) = 4Δf Σ a[k,:]' inv(S[k,:,:]) b[k,:]`.
  The first argument is conjugated, as with Julia's `dot`.
- The real inner product is `real(C(a,b))`; optimal SNR squared is `(h|h)`.
- The complex matched-filter scalar is `C(h,d)/sqrt((h|h))`. A zero template
  has zero optimal SNR but raises `DomainError` for matched-filter SNR.
- `residual_squared(d,h,m)` computes `(d-h|d-h)`.
- `log_likelihood_ratio(h,d,m)` computes `real(C(h,d)) - (h|h)/2` for a fixed
  PSD, relative to the noise-only hypothesis. It is not a normalized density.
- Uniform increasing nonnegative `AbstractRange` grids only. No interpolation,
  resampling, time shifting or channel reordering occurs automatically.
- Frequency bounds are inclusive: select bins `f >= minimum_frequency` and
  `f <= maximum_frequency`. DC is always excluded. `nsamples` identifies and
  excludes even-length Nyquist on a complete real-FFT grid; odd-length final
  bins remain active. Without it the last positive bin is treated normally.
- Active PSDs must be finite and strictly positive. Out-of-band zeros are
  ignored before division/factorization. When using `psd(::GaussianNoise)`, set
  the same analysis band and pass `nsamples` to exclude its zero endpoints.

All channels in a joint metric share a grid. Independent groups with different
grids can be evaluated separately and their contractions added. Correlated
channels must first be represented on a common Fourier basis. When dropping a
channel, construct the retained covariance submatrix and refactor it; slicing
the inverse covariance would instead describe a conditional noise model.

## Differentiation and storage

Validation and factorization run once on the host. Evaluation preserves strain
element types, including ForwardDiff duals. The intended initial AD contract is
derivatives through strain/waveform parameters with a fixed metric; differentiating
PSD/CSD preparation is outside scope. For compiled matched-filter SNR use
`matched_filter_snr(h, d, metric; check=false)` after establishing a nonzero
template on the host; the default throws for zero templates and has a
value-dependent validation branch. Heavy backend acceptance tests are in
`integration/GWInference`.

CPU values and gradients have been checked with Reactant 0.2.295, Enzyme
0.13.216 and Mooncake 0.5.54 on Julia 1.12.0, including waveform → projection →
correlated likelihood composition. GPU execution has not been tested. The
integration example materializes the fixed projection frequency vector once to
avoid differentiating Julia's range-broadcast machinery.

Strain arrays use `(frequency, channel)` storage; internal whitening blocks use
`(output_channel, input_channel, frequency)`. The [layout benchmark](../../benchmarks/GWInference)
records the measured decision. Native host contractions use allocation-free
StaticArrays kernels; traced arrays use frequency-wise broadcast/reductions.
Both flattening orders are mathematically equivalent when the metric uses the
same permutation; evaluation never materializes the global flattened matrix.

## Development

```sh
julia --project=packages/GWInference -e 'using Pkg; Pkg.test()'
just test-integration-inference
just benchmark-inference
```

Tests compare against an explicit dense metric, analytic correlated models,
Gaussian-noise normalization, ForwardDiff gradients and optional bilby fixtures.
Run `uv run --project python python python/generate_inference.py` for fixtures;
CI requires them. Heavy integration tests additionally exercise Enzyme, Mooncake,
Reactant CPU and waveform → detector projection → likelihood gradients.
