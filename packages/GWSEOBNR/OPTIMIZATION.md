# ForwardDiff optimization plan

Status: **baseline and profile measured; candidates 1 and 2 rejected, nothing applied.** ForwardDiff is the supported
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

## Profile

Sampling profile of the 7-parameter gradient at M=40, q=1, chunk 7, 5 gradient calls
(3143 samples). Percentages are the share of samples with the function anywhere on the
stack, so they overlap:

| function | share |
| --- | --- |
| `sample` (per-output-sample jets, `fd.jl`) | 63% |
| `eob_rhs` (called from the ODE solve and from dense output) | 60% |
| `radiation_flux` | 39% |
| `inspiral_mode` | 38% |
| `evolve` (ODE solve, including its step stages) | 36% |
| `time_jet` | 29% |
| `ode_interpolation` (dense output) | 28% |
| `peak_time` | 1% |
| `real_fft` | 0.2% |
| `initial_conditions` | 0.1% |

The cost is in the per-sample evaluation, and within it in dense-output interpolation,
which for Vern9 lazily computes extra stages with further `eob_rhs` calls on Duals.
`eob_rhs`, and `radiation_flux` inside it, is the hot kernel for both the solve and the
samples.

## Candidates, in order of expected value

1. **Redundant flux coefficients.** Tried and rejected. Building `c` and the flux
   coefficients once per draw in `initial_solution`, instead of on every Newton residual
   evaluation, changed chunk-7 allocations by under 0.2% (1585 to 1583 MiB, 688 to
   690 MiB) and the timing difference was inside run-to-run noise. The profile below shows
   why: the whole initial-condition solve is 0.1% of the samples.
2. **FFT adapter.** Not worth doing. `real_fft` is 0.2% of the samples in the profile
   below, so batching the transforms cannot gain anything measurable.
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
