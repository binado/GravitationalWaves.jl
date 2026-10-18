# Frequency-domain metric layout benchmark

Run from the repository root with Julia 1.12 or newer:

```sh
julia --project=benchmarks/GWMatchedFiltering -e 'using Pkg; Pkg.instantiate()'
julia --threads=1 --project=benchmarks/GWMatchedFiltering benchmarks/GWMatchedFiltering/layouts.jl > /tmp/gw-inference-layouts.csv
```

The script checks numerical equivalence before measuring. It reports median
latency, allocated bytes, and allocation counts after warmup, with one BLAS
thread. Julia version and CPU model appear in the output. Construction is
measured separately from repeated evaluation. The pipeline measurements include
waveform generation, projection, and any required transposition.

Compare frequency-first and channel-first storage, plus frequency-first strains
with contiguous whitening blocks, using scalar and StaticArrays block kernels,
diagonal contractions, and the production kernels. Test
1, 3, and 6 channels at 1,025, 8,193, and 32,769 active frequencies. Treat the
small-kernel timings as diagnostic; choose storage using representative repeated
statistics and pipeline costs, and check the AD integrations before adopting a
different kernel.

## Storage decision

Use `(frequency, channel)` strains and `(output_channel, input_channel, frequency)`
internal whitening factors. This keeps projected detector spectra contiguous
without a transpose, while each small whitening matrix is contiguous. Public
spectral-density input remains `(frequency, channel, channel)`; preparation
reorders it once. Native host kernels use StaticArrays and allocate no memory for
inner products, norms, or matched-filter scalars. Traced arrays use reductions
over frequencies instead of scalar unrolling.

Measured with Julia 1.12.0 on a Sapphire Rapids CPU, one Julia/BLAS thread.
Representative median static-block contraction times (microseconds), with zero
allocations in each case:

| Channels | Frequencies | Frequency-first throughout | Channel-first throughout | Frequency-first strain, contiguous blocks |
| --- | --- | --- | --- | --- |
| 3 | 8,193 | 146 | 152 | 142 |
| 3 | 32,769 | 609 | 634 | 596 |
| 6 | 8,193 | 558 | 481 | 489 |
| 6 | 32,769 | 2,288 | 1,949 | 2,017 |

The mixed storage wins the three-channel case without input conversion. Full
channel-first storage can be faster for six channels; the chosen convention
balances the main three-detector workload with existing projection output. The
three-channel waveform/projection benchmark also showed no benefit from
transposing strains (approximately 0.84 ms versus 0.88 ms in the final run).

See [raw comparison results](layout_comparison.csv) and the
[full benchmark snapshot](results.csv). The final run had no other Julia jobs
running. Contiguous blocks improved this static kernel by about 2–13% over
frequency-first blocks; an earlier run with background compilation overstated
that difference. These are host-specific
measurements, not CI timing thresholds. Rerun on deployment hardware; compiled
Reactant performance is a separate backend measurement. The initial allocating
broadcast implementation was substantially slower on the host, motivating the
native specialization and the StaticArrays runtime dependency.

## GWSeries delegation

The layout benchmark also records `series/blocks` and `series/diagonal` alongside
production array calls. For a focused comparison, including native matrix views:

```sh
julia --threads=1 --project=benchmarks/GWMatchedFiltering benchmarks/GWMatchedFiltering/series.jl
```

This validates equal numerical results before timing, with fixed epoch/channel
metadata and three channels at 8,193 active frequencies. Metadata validation is
included in series timings; metric preparation and series construction are
excluded. FFTs, storage adaptation, and whole-series compilation are separate
contracts. Rerun on deployment hardware; these measurements are not CI thresholds.

Measured on Apple M4 with Julia 1.12.7, one Julia/BLAS thread and no competing
Julia jobs. Median microseconds, with zero bytes and zero allocations throughout:

| Metric | Storage | Arrays | Series |
| --- | --- | --- | --- |
| Diagonal | Matrix | 15.000 | 15.042 |
| Diagonal | Matrix view | 15.041 | 15.042 |
| Correlated blocks | Matrix | 26.417 | 26.459 |
| Correlated blocks | Matrix view | 26.375 | 26.417 |

Observed wrapper overhead was 0–0.3% in these cases. See
[raw series results](series_results.csv). This establishes representative host
behavior, not a guarantee for other storage, sizes, hardware, or compiled backends.
