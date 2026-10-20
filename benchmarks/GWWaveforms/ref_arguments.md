# Ref arguments versus anonymous functions

```sh
julia --threads=1 --project=integration/GWWaveforms benchmarks/GWWaveforms/ref_arguments.jl
```

Instantiate both the integration and benchmark environments first, as described in
[the gradient benchmark](gradients.md). The script retains the previous anonymous-function
assembler and compares it directly with the public production API, which now broadcasts the
scalar kernels with a shared `Ref(model)`, `Ref(scale)` and `Ref(cutoff)`.

All implementations prepare the model once and return separate amplitude/phase vectors.
Outputs are checked for exact equality before timing. Measurements use Julia 1.12.7, Apple M4,
one Julia thread, the same Float64 parameters as the earlier benchmarks, and a 20–1000 Hz
configured waveform interval on a grid spanning 20–1024 Hz. Tests and formatting processes
finished before these benchmarks began. Compilation and AD cache preparation are excluded.

Primal warmed medians over 200 single evaluations:

| Model | Bins | Anonymous (μs) | Ref (μs) | Anonymous bytes / allocations | Ref bytes / allocations |
|---|---:|---:|---:|---:|---:|
| IMRPhenomD | 4,017 | 97.29 | 137.00 | 66,432 / 12 | 65,664 / 6 |
| IMRPhenomD | 32,129 | 770.17 | 1,101.83 | 525,184 / 12 | 524,416 / 6 |
| IMRPhenomXAS | 4,017 | 187.50 | 225.33 | 66,432 / 12 | 65,664 / 6 |
| IMRPhenomXAS | 32,129 | 1,471.71 | 1,792.40 | 525,184 / 12 | 524,416 / 6 |

Prepared gradients of the same seven-parameter, phase-sensitive loss as the earlier gradient
benchmark, at 4,017 frequencies. These use equivalent Float64-backed ranges to avoid the
Mooncake metadata-copy limitation described there. Timings are medians of 100 single calls:

| Model | Backend | Anonymous (μs) | Ref (μs) | Anonymous bytes / allocations | Ref bytes / allocations |
|---|---|---:|---:|---:|---:|
| IMRPhenomD | ForwardDiff | 646.73 | 998.29 | 2,621,824 / 18 | 2,621,824 / 18 |
| IMRPhenomD | Mooncake | 2,332.73 | 2,616.60 | 671,152 / 328 | 664,544 / 231 |
| IMRPhenomXAS | ForwardDiff | 859.56 | 1,326.42 | 4,934,048 / 9,906 | 4,929,312 / 9,900 |
| IMRPhenomXAS | Mooncake | 2,281.69 | 2,504.56 | 669,408 / 229 | 662,672 / 130 |

In this comparison, Ref lowers the primal allocation count from 12 to 6, saving a fixed
768 bytes per call. Primal timings increase by about 41% for PhenomD and 20% for PhenomXAS at
4,017 frequencies. ForwardDiff timings increase by about 54% for both models; Mooncake timings
increase by about 12% and 10%, respectively. The counts improve more than the allocated bytes.
These measurements describe these implementations and do not establish a general performance
rule about Ref; profiling would be needed to attribute the slowdown.

All 54 gradient checks passed, including finite differences and cache reuse across equal
masses, swapped components and a large cutoff movement. ForwardDiff gradients agree exactly
with the scalar-loop reference at the benchmark points; Mooncake's largest relative L2 error
is `6.294e-15`, unchanged by the Ref conversion. The production Ref change also passes all
22,863 package tests and 134 AD integration checks.

The complete machine-readable results are [primal](ref_primal_results.csv) and
[gradients](ref_gradient_results.csv). These are same-session comparisons; timing ratios
should be calculated within this report rather than across separate benchmark runs.
