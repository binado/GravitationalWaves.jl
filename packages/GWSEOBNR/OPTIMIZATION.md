# ForwardDiff optimization plan

Status: **baseline measured; no optimization applied yet.** ForwardDiff is the supported
backend, and the Mooncake extension reuses it: its `spectrum_vector` rule runs one
ForwardDiff Jacobian pass (`GWSEOBNRMooncakeExt`). A faster Dual path therefore speeds
up both.

## Baseline

`just benchmark-seobnr` (`benchmarks/GWSEOBNR/gradient.jl`) measures the primal and the
7-parameter ForwardDiff gradient of `evaluated_spectrum` on the `0:0.125:1024` grid
(8193 bins) at chunk sizes 1, 2, 4 and 7. Julia 1.12.7, Apple M4, one thread, minimum of
5 runs after a warm-up call:

| point | primal | chunk 1 | chunk 2 | chunk 4 | chunk 7 |
| --- | --- | --- | --- | --- | --- |
| M=40, q=1 | 0.33 s / 356 MiB | 3.18 s / 3907 MiB | 2.10 s / 2725 MiB | 1.04 s / 2020 MiB | 0.80 s / 1585 MiB |
| M=60, q=2, spinning | 0.12 s / 151 MiB | 1.09 s / 1686 MiB | 0.75 s / 1176 MiB | 0.50 s / 875 MiB | 0.40 s / 688 MiB |

Allocations repeat exactly between runs. Wall time does not: the M=40 primal ranged from
0.24 s to 0.33 s over two runs on the same machine, so compare changes by allocations
first and re-run before trusting a timing difference under about 20%.

Chunk 7 is the fastest and is also ForwardDiff's default for 7 inputs, so candidate 5
below needs no change.

## Rules

- Run `just benchmark-seobnr` (with `JULIA="julia +1.12"` to pin the version) before and
  after each change and compare against the table above.
- Keep a change only if it measurably wins and leaves the acceptance gates unchanged
  (130 passing and 13 failing on Julia 1.12.7, the same gates failing; see
  `ACCEPTANCE.md`). Run `just test-integration-seobnr` and
  `just test-integration-seobnr-mooncake` before and after.
- Do not change the physics. The remnant fits, flux and Hamiltonian follow
  pySEOBNR operation order, and the fidelity gates depend on it.
- Measure on one Julia version. Results differ between 1.12 and 1.13, so state which.

## Candidates, in order of expected value

1. **Redundant flux coefficients.** `flux_coefficients` fills eight `zeros(T, 9, 9, 11)`
   arrays in the Dual element type. `initial_residual` (`src/initial.jl`) calls it on
   every Newton residual evaluation, and `check_initial_solution` calls it again through
   the Jacobian, while `evolve` computes it once. Compute it once per draw and pass it in.
   The coefficients depend on the differentiated parameters (`nu`, `chi_S`, `chi_A`),
   so they cannot be computed in `Float64`. Options are reuse, preallocation, or
   `SArray` storage if the sizes allow it.
2. **FFT adapter.** `real_fft` in `src/fft.jl` makes one `rfft` per partial, `P + 1`
   transforms per chunk. A single batched `rfft` over a matrix of the value and
   partials, with a preplanned FFT, should do the same work with less overhead.
3. **Per-sample jets.** `sample(j)` in `src/fd.jl` evaluates five-point `time_jet`s and
   mode evaluations for every output sample. Share the Dual evaluation across samples
   where the inputs coincide, and avoid growing `candidates` in `peak_time`.
4. **ODE error control.** Check whether the default `internalnorm` includes the
   partials of the Dual state. If so, step rejections differ from the `Float64`
   sequence. A norm on primal values only would match the `Float64` step sequence, but
   it must not loosen derivative accuracy; the derivative convergence gate
   (`derivative_convergence`) decides this.
5. **Chunk size.** Resolved by the baseline: the default is already 7, which is the
   fastest. The Mooncake rule has the same 7 inputs; confirm it uses that default.
6. **Type stability.** Run `@code_warntype` and JET on the Dual path, in particular
   closures that capture `sol` and the `Float64(primal(...))` casts in `src/roots.jl`.

## Mooncake-specific notes

- The rule's cost is one Jacobian pass plus a real matrix of `4 nf × 7` entries held for
  the pullback (`nf` frequency bins, real and imaginary parts of both polarizations). On
  the `0:0.0625:2048` acceptance grid, `nf = 32769`, that is about 7 MB. Check the
  memory if grids grow.
- Compile time for differentiating the flux code directly with Mooncake was about 35
  minutes (`diagnostics/mooncake-2026-10-19.toml`). Items 1 and 3 may reduce it, but the
  rule makes it irrelevant for `hphc_fd` and likelihood gradients.

## Out of scope

Fixing the NQC and spectral-L² fidelity gates so `is_differentiable` can become true is
a separate task.
