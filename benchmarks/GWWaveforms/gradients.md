# Gradient comparison

```sh
julia --project=benchmarks/GWWaveforms -e 'using Pkg; Pkg.instantiate()'
julia --project=integration/GWWaveforms -e 'using Pkg; Pkg.instantiate()'
julia --threads=1 --project=integration/GWWaveforms benchmarks/GWWaveforms/gradients.jl
```

[The complete results](gradient_results.csv) include loops, tuple broadcasts and separate
broadcasts at 503 and 4,017 frequencies. The script loads BenchmarkTools from the benchmark
environment and the AD packages from the integration environment; runtime dependencies are unchanged.

Measurements use Julia 1.12.7, Apple M4, one Julia thread, ForwardDiff 1.4.6 and Mooncake 0.5.64.
The masses, spins, distance, inclination and phase are `(36, 29, 0.5, 0.3, 500, 0.4, 0.2)`.
Frequencies span 20–1024 Hz, with a configured waveform interval of 20–1000 Hz and the usual
intrinsic cutoffs. The table below uses spacing 0.25 Hz (4,017 frequencies).

The differentiated scalar loss depends on both polarizations and their phases:

```julia
h = (1e22 .* amplitude) .* cis.(phase)
plus = ((1 + cos(inclination)^2) / 2) .* h
cross = (-im * cos(inclination)) .* h
loss = sum(abs2, plus) + real(sum(cross .* conj.(fixed_template)))
```

The fixed template uses the same parameters with `mass_1` increased by 1%. Each implementation
performs identical parameter conversion, coefficient preparation and downstream loss evaluation.
The seven gradient coordinates are component masses, two aligned spins, distance, inclination
and reference phase. Sampling and generator configuration are fixed.

ForwardDiff uses `gradient!` with a preallocated result/configuration and chunk size 7.
Mooncake uses `value_and_gradient!!` with a prepared cache. Times are warmed medians of 100
single evaluations. Allocation measurements cover each call, excluding compilation, configuration
construction and Mooncake cache preparation/retained memory. Both methods evaluate the primal
loss while computing its gradient; Mooncake also returns that loss value.

| Model | Backend | Scalar loop (μs) | Separate broadcasts (μs) | Change |
|---|---|---:|---:|---:|
| IMRPhenomD | ForwardDiff | 612.88 | 577.25 | 5.8% faster |
| IMRPhenomD | Mooncake | 1,938.13 | 2,274.88 | 17.4% slower |
| IMRPhenomXAS | ForwardDiff | 883.71 | 900.48 | 1.9% slower |
| IMRPhenomXAS | Mooncake | 2,075.83 | 2,247.69 | 8.3% slower |

Per-call allocated bytes and allocation counts:

| Model | Backend | Scalar loop bytes / count | Separate broadcasts bytes / count |
|---|---|---:|---:|
| IMRPhenomD | ForwardDiff | 2,621,824 / 18 | 2,621,824 / 18 |
| IMRPhenomD | Mooncake | 663,888 / 228 | 671,152 / 328 |
| IMRPhenomXAS | ForwardDiff | 4,929,312 / 9,900 | 4,934,048 / 9,906 |
| IMRPhenomXAS | Mooncake | 662,016 / 127 | 669,408 / 229 |

All 154 checks passed. ForwardDiff gradients agree exactly across implementations at the
benchmark points; the largest Mooncake relative L2 gradient error versus the loop's
ForwardDiff result is `6.294e-15`. Five-point finite differences independently check the
coarse-grid baseline (`rtol = 1e-5`). Implementation gradients agree at additional points
with swapped components, equal masses and masses `(180, 120)` moving many cutoff bins.
The same Mooncake cache is reused across these points, with agreement against ForwardDiff
checked at `rtol = 1e-9`. These checks evaluate derivatives at the chosen points and do not
assert differentiability across jumps of the hard frequency mask.

The ordinary `20.0:df:1024.0` range triggered a Mooncake 0.5.64 cache-preparation error while
copying captured `Base.TwicePrecision` metadata on a subsequent broadcast cache. This benchmark
uses `StepRangeLen(20.0, df, n)` with Float64 metadata instead. For the binary spacings tested,
the script verifies exact equality of all sampled frequencies against the ordinary range.
No custom derivative rules or changes to Mooncake are needed. The earlier primal timing
results use ordinary ranges; the gradient table compares all variants using the same
Float64-backed range.

Gradient performance does not follow the primal timing ratios: separate broadcasts modestly
help PhenomD with ForwardDiff here, while both Mooncake comparisons favor the loop. These
measurements do not identify the compiler or AD internals responsible for each difference.
