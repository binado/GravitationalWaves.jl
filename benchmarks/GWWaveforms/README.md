# Scalar Phenom kernels

The separate-broadcast measurements below used anonymous functions. The later
[`Ref` comparison](ref_arguments.md) measures the current production implementation
against that version.


```sh
julia --project=benchmarks/GWWaveforms -e 'using Pkg; Pkg.instantiate()'
julia --threads=1 --project=benchmarks/GWWaveforms benchmarks/GWWaveforms/scalar_kernels.jl
```

The benchmark compares four implementations:

- The direct-frequency loops before scalar-kernel extraction, saved in `previous_loops.jl`.
- The previous scalar loops calling scalar amplitude/phase kernels.
- Broadcasting tuple-returning scalar kernels and extracting their components into output arrays.
- Broadcasting scalar amplitude and phase kernels separately into the two output arrays,
  without an intermediate array of tuples. Phase canonicalization consumes the signed
  amplitudes before their absolute values are taken in place.

All variants prepare coefficients once per waveform, return separate amplitude and phase
vectors, and must produce exactly equal outputs before timing. The loop and tuple-broadcast
alternatives live in the benchmark script; separate broadcasts now call the public production
method directly.

Timings include parameter conversion, model preparation and output allocation, excluding
compilation and host validation. Each specialization is warmed up before 200 samples with
one evaluation per sample. The script reports median/minimum times, bytes and allocations.

Measured on 2026-10-20 with Julia 1.12.7, Apple M4, one Julia thread, Float64 parameters
(masses 36 and 29 solar masses, spins 0.5 and 0.3), and uniform ranges from 20 to 1024 Hz.
The configured waveform interval is 20–1000 Hz; intrinsic model cutoffs still apply.

Median times in microseconds:

| Model | Bins | Previous loop | Scalar loop | Tuple broadcast | Separate broadcast |
|---|---:|---:|---:|---:|---:|
| IMRPhenomD | 503 | 9.46 | 9.67 | 11.77 | 12.67 |
| IMRPhenomD | 4,017 | 73.63 | 71.54 | 88.96 | 94.21 |
| IMRPhenomD | 32,129 | 580.25 | 566.52 | 679.42 | 747.96 |
| IMRPhenomXAS | 503 | 24.44 | 24.42 | 27.13 | 25.08 |
| IMRPhenomXAS | 4,017 | 179.46 | 178.96 | 197.63 | 182.52 |
| IMRPhenomXAS | 32,129 | 1,413.92 | 1,411.88 | 1,562.60 | 1,445.60 |

Allocated bytes / allocation counts (the same for both models):

| Bins | Previous and scalar loops | Tuple broadcast | Separate broadcast |
|---:|---:|---:|---:|
| 503 | 8,320 / 6 | 17,328 / 14 | 9,072 / 11 |
| 4,017 | 65,664 / 6 | 132,032 / 15 | 66,432 / 12 |
| 32,129 | 524,416 / 6 | 1,033,152 / 15 | 525,184 / 12 |

Separate broadcasts eliminate the grid-sized temporary; the extra memory is now small and
fixed rather than proportional to grid size. Broadcasting does not require a tuple array.
The original tuple-broadcast measurements therefore do not establish that broadcasting
in general doubles memory usage.

PhenomD's separate amplitude and phase passes each compute `useful_powers(Mf)`, whereas the
combined scalar kernel shares that calculation. Separate broadcasting is about 32% slower
than the scalar loop at 4,017 and 32,129 bins here. PhenomXAS shares less frequency-dependent
work and its separate broadcasts are within about 2–3% of the loop. These are measurements
of complete implementations, not a profiling-based attribution of every timing difference.
Production uses separate broadcasts for the preferred Julia style, accepting the measured
timing and allocation tradeoffs. The script retains both loop baselines for comparison.

The [gradient comparison](gradients.md) measures ForwardDiff and Mooncake timings,
allocations, numerical agreement and reverse-cache reuse for these same implementations.
