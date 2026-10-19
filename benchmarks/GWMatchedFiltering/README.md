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

## Domain-first inference

For a focused comparison of domain-first evaluation against the underlying array
contraction kernel, including native matrix views:

```sh
julia --threads=1 --project=benchmarks/GWMatchedFiltering benchmarks/GWMatchedFiltering/domains.jl
```

The comparison uses one and three channels, vector/matrix storage and their views,
and 8,193 active frequency bins. It validates data
and metric compatibility before timing, and checks equal numerical results. Both
paths include the same physical `4df` normalization. No compatibility checks or
coordinate-array construction occur during evaluation. Metric/domain construction
is excluded. Measurements are hardware-specific, not CI timing thresholds.

Measured after integration checks on Apple M4, Julia 1.12.7, one Julia/BLAS
thread. Median microseconds; every case allocated zero bytes and zero objects:

| Channels | Metric | Storage | Kernel | Domain-first |
| --- | --- | --- | --- | --- |
| 1 | Diagonal | Vector | 5.000 | 5.000 |
| 1 | Diagonal | Vector view | 5.000 | 5.000 |
| 1 | Blocks | Vector | 5.584 | 5.584 |
| 1 | Blocks | Vector view | 5.625 | 5.625 |
| 3 | Diagonal | Matrix | 15.000 | 15.000 |
| 3 | Diagonal | Matrix view | 15.000 | 15.041 |
| 3 | Blocks | Matrix | 26.417 | 26.375 |
| 3 | Blocks | Matrix view | 26.375 | 26.375 |

Domain-first timings were within 0.3% of the raw kernel in these measurements,
including small variations in both directions. See [raw results](domains_results.csv).
This is representative host evidence, not a guarantee for other sizes, storage,
hardware, or compiled backends. Previous series-wrapper measurements have been retired.
