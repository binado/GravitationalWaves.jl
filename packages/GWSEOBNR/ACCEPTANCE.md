# Acceptance pending

The prototype is **not accepted** against the implementation plan. Keep
`is_differentiable(SEOBNRv5HM()) == false`. The numerical and AD gates are
unchanged, and all nine cases remain in the acceptance runner.

Normal package tests pass on Julia 1.12.7: Hamiltonian values and analytic
first/second derivatives, full multipolar flux, calibration and remnant fits,
initial conditions, raw (2,2) modes, QNM interpolants for both spin signs, local
FFT/tail AD adapters, root failures, coverage/attachment failures, output masking,
and a compact short-waveform oracle. The CI-selection unit tests (15 checks),
Runic across every Julia file, and all four existing heavy AD suites (1,350
checks) also pass.

The complete SEOBNR acceptance run passes intrinsic directional derivatives at
two adjacent perturbations, extrinsic finite differences, waveform/derivative
solver convergence, and the offset H1/L1 fixed-PSD likelihood gradient. It still
fails the NQC coefficient fidelity gate in all nine cases and the primary
spectral L² gate in four cases. Raw, unoptimized mismatches pass in every case.

The Julia 1.12.7 run after the production NonlinearSolve migration reports
**130 passing and 13 failing checks**, with no test errors. The committed
[diagnostic record](diagnostics/acceptance-production-nonlinearsolve-2026-10-19.toml)
contains solver steps, timing roots, matching
coefficients, raw/optimized errors, five-point finite-difference sweeps,
convergence, and allocation/timing measurements. Large raw and conditioned
Python fixtures remain under gitignored `reference/seobnr`.

| Total mass (M☉) | q | Relative spectral L² | Unoptimized mismatch | Relative NQC error |
|---|---|---|---|---|
| 40 | 1 | 1.478e-04 | 1.464e-08 | 3.336e-03 |
| 40 | 2 | 8.098e-04 | 3.545e-07 | 7.469e-02 |
| 40 | 4 | 6.865e-05 | 3.362e-09 | 2.569e-04 |
| 60 | 1 | 4.746e-04 | 1.608e-07 | 7.154e-03 |
| 60 | 2 | 1.095e-03 | 9.162e-07 | 7.448e-02 |
| 60 | 4 | 2.481e-04 | 6.195e-08 | 2.834e-04 |
| 100 | 1 | 1.733e-03 | 7.746e-07 | 7.222e-03 |
| 100 | 2 | 1.230e-03 | 6.437e-07 | 3.839e-02 |
| 100 | 4 | 1.001e-03 | 1.767e-07 | 9.445e-05 |

Required limits are respectively **1e-3**, **1e-5**, and **1e-9**. Diagnostic
optimization of time/phase does not change any primary comparison.

The native NQC implementation currently differentiates continuous Vern9 dense
output using a fixed five-point time stencil. The pinned upstream construction
uses local sampled spline derivatives. These matching differences require
further work; they are not covered by a claim that only physical timing changed.
The spectrum failures must be resolved against the pinned oracle before model
acceptance or enabling the differentiability capability flag.

Reproduce with `just fixtures-seobnr` followed by
`just test-integration-seobnr`. The latter intentionally exits unsuccessfully
while any gate remains unmet and writes diagnostics even when tests fail.

The initial-condition migration leaves the identities of all 13 failed gates
unchanged. The largest change in the spectral-L² error statistic versus the
[pre-migration record](diagnostics/acceptance-2026-10-19.toml) is `2.64e-8`.
Waveform and derivative convergence changes remain below `2.84e-8` and
`4.08e-6`, respectively; the offset likelihood-gradient finite-difference error
is `6.58e-8`. These regression results do not resolve the existing fidelity gates.
