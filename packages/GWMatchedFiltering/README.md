# GWMatchedFiltering

Joint frequency-domain noise metrics, inner products, and SNRs for ordinary
arrays. Requires GWTimeFrequencyDomains, LinearAlgebra, and StaticArrays.

```julia
using GWTimeFrequencyDomains, GWMatchedFiltering

fd = FrequencyDomain(; frequencies=0:0.25:1024, epoch=1126259460.0,
                     channels=(:H1, :L1))
spectrum = fill(1e-46, nsamples(fd), 2)
metric = FrequencyDomainMetric(fd, spectrum; minimum_frequency=20, nsamples=8192)
# h and data are caller-supplied (frequency, channel) arrays.
check_data(fd, h)
check_data(fd, data)
check_compatible(fd, metric)
check_template(fd, h, metric)

# Repeated inference evaluation: no metadata or template-value validation.
optimal_snr(fd, h, metric)
matched_filter_snr(fd, h, data, metric)
noise_weighted_inner_product(fd, h, data, metric)
likelihood_statistics(fd, h, data, metric)
log_likelihood_ratio(fd, h, data, metric)
```

Arrays use `(frequency, channel)` layout; a vector denotes one channel. Keep
independently sourced operands' domains and call `check_same_grid` before combining
them. Passing one domain to evaluation asserts that both arrays use its grid,
epoch, and channel order. Detached arrays cannot verify that assertion.

## Prepared noise metrics

`FrequencyDomainMetric(domain, spectrum; minimum_frequency, maximum_frequency,
nsamples=nothing)` validates and prepares numerical state on the host. PSDs use
shape `(frequency, channel)` or a vector for one channel. Correlated one-sided
spectral matrices use `(frequency, channel, channel)`, with
`S[k,i,j] = 2df * E[n[k,i] * conj(n[k,j])]` for interior bins.

The prepared metric records its frequency grid/channel order, active bins, and
weights/whitening factors. It has no epoch requirement and stores no redundant
normalization: evaluation takes `4frequency_spacing(domain)` from the domain.
`check_compatible(domain, metric)` checks its setup grid and channel order.

Band bounds are inclusive. DC is always excluded. Supply `nsamples` on a complete
zero-origin FFT grid to exclude even-length Nyquist; odd-length final bins remain
active. Without it the last positive bin is ordinary. Active PSDs must be finite
and positive. Active spectral matrices must be finite, exactly Hermitian, and
positive definite. Out-of-band zeros are allowed; no interpolation, jitter,
pseudoinverse, resampling, or channel reordering is performed.

Correlated preparation caches inverse lower Cholesky factors `W[k] = L[k] \\ I`.
Contractions use `(W*a)' * (W*b)` without constructing a global covariance or
inverse. For channel subsets, prepare the retained covariance submatrix again;
slicing inverse covariance would describe a conditional model.

## Statistics and inference preconditions

All evaluation methods take a `FrequencyDomain` first:

- `complex_inner_product(fd, a, b, metric)` computes `4df Σ a' inv(S) b`.
- `noise_weighted_inner_product` takes its real part.
- `optimal_snr_squared` and `optimal_snr` compute `(h|h)` and its square root.
- `matched_filter_snr` computes the complex scalar `C(h,data)/sqrt((h|h))` at
  the supplied alignment; no time/phase maximization is performed.
- `residual_squared(fd, data, h, metric)` computes `(data-h|data-h)`.
- `likelihood_statistics` returns `(h_inner_d, h_inner_h)` contractions.
- `log_likelihood_ratio` computes `real(C(h,data)) - (h|h)/2` for fixed PSDs.

`check_template` checks domain/metric compatibility, array shape, and finite
positive template norm, throwing for zero/nonfinite templates. Evaluation does
not invoke it: callers must enforce changing parameter support and nonzero norms.
Zero templates passed directly to matched-filter evaluation produce nonfinite
results. Shape/domain mismatch is a caller error after validation, not a checked
operation in an inference loop.

## Differentiation and storage

Fixed domains and prepared metrics are captured during compilation. Runtime
samples and physical parameters are differentiated; metric preparation and
metadata are outside the AD contract. Metric kernels and the separately validated
inference boundary are checked with ForwardDiff, Enzyme, Mooncake, and Reactant CPU.
Waveform → domain-first projection → correlated likelihood composition is checked
with ForwardDiff and Mooncake, matching GWWaveforms' supported AD backends.
Dynamic whole-domain Reactant arguments and GPU execution have not been established.

Native host matrices and their views retain the allocation-free contraction
specializations. Traced arrays use frequency-wise broadcast/reduction kernels.
See the [benchmarks](../../benchmarks/GWMatchedFiltering) for measured comparisons.

## Migration to 0.2

`FrequencyDomainStrain`, `TimeSeries`, and `FrequencySeries` are not part of this
API. Keep arrays unchanged and construct a detached `FrequencyDomain`. Replace
`statistic(h, data, metric)` with `statistic(fd, h, data, metric)`, construct metrics
from `(fd, spectrum)`, and run checks before inference. The `matched_filter_snr`
`check` keyword has been removed; use `check_template` explicitly.

Run `julia --project=packages/GWMatchedFiltering -e 'using Pkg; Pkg.test()'`,
`just test-integration-inference`, and `just benchmark-inference`. Generate bilby
fixtures with `uv run --project python python python/generate_inference.py`;
CI requires them. Formatting uses Runic.
