# GWMatchedFiltering

Joint frequency-domain inner products, SNRs and fixed-PSD Gaussian likelihoods.
Series use `(frequency, channel)` storage and one-sided physical Fourier normalization.

```julia
using GWSeries, GWMatchedFiltering
f = 0:0.25:1024
h = FrequencySeries(ones(ComplexF64, length(f), 2); frequency_spacing=step(f),
                    epoch=1126259460.0, channels=(:H1, :L1))
data = copy(h)
metric = FrequencyDomainMetric(h, ones(length(f), 2);
                             minimum_frequency=20, nsamples=8192)
log_likelihood_ratio(h, data, metric)
```

Prepare without a prototype using `FrequencyDomainMetric(spectrum; frequencies,
channels, minimum_frequency, maximum_frequency, nsamples=nothing)`.
Independent PSDs have shape `(frequency, channel)` (or a vector for one channel).
Correlated spectra have shape `(frequency, channel, channel)`, with
`S[k,i,j] = 2Δf E[n[k,i] conj(n[k,j])]`. Active spectral blocks must be exactly
Hermitian and positive definite. Preparation caches whitening factors, never an
inverse covariance matrix. Diagonal PSDs must be positive and finite in the active band.

Band limits are inclusive. DC is excluded. Supplying the original time sample
count also excludes even-length Nyquist; odd-length last bins remain active.
Stationary metric compatibility checks grid and ordered channels but ignores epoch.
FrequencySeries always describes a complete one-sided grid. No interpolation or channel reordering is implicit.

Series methods check alignment and metric compatibility before unwrapping storage.
Raw-array methods with the same argument order perform no metadata checks:

```julia
check_same_grid(h, data)
check_compatible(h, metric)
check_template(h, metric)       # finite positive norm, required for matched-filter SNR
ha, meta = flatten(h)
da, _ = flatten(data)
log_likelihood_ratio(ha, da, metric) # suitable for repeated compiled work
```

Available statistics are `complex_inner_product(a,b,metric)`,
`noise_weighted_inner_product`, `matched_filter_snr`, `residual_squared`,
`likelihood_statistics`, `log_likelihood_ratio`, `optimal_snr_squared(h,metric)`,
and `optimal_snr(h,metric)`. The inner product conjugates its first argument and
uses `4Δf` normalization supplied by the prepared metric.

Native matrices/vectors and their views retain allocation-free contraction
kernels; traced arrays retain broadcast/reduction kernels. Heavy AD acceptance
lives in `integration/GWMatchedFiltering`. Construct metrics outside differentiation;
only strain parameters are active. Single-channel raw vector kernels avoid
reshaping mixed-activity array headers inside Enzyme calls.

For band-only work, keep `bins = frequency_bins(h; fmin=20, fmax=512)` and prepare
`FrequencyDomainMetric(psd_values[bins, :]; frequencies=frequencies(h)[bins], channels=channels(h))`.
Pass ordinary arrays on that same band to its statistics. The series constructor does not
represent cropped grids. Explicit `nsamples` is only accepted with a complete zero-origin grid.
