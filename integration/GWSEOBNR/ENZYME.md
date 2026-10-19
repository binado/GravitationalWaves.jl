# Exploratory Enzyme reverse-mode probes

Measured on 2026-10-19, macOS arm64, Julia 1.12.7, Enzyme 0.13.216.
Native reverse gradients work for several physics components. Complete waveform
and detector likelihood reverse-mode support remains unavailable. These are
exploratory probes at one physical draw, not the nine-case acceptance suite.

The draw is `(mass_1, mass_2, chi_1, chi_2, distance, inclination, phase)` =
`(40, 20, 0.4, -0.3, 500, 0.7, 0.3)`, with masses in solar masses and distance
in Mpc. Complete spectra use `0:0.125:1024` Hz and default solver tolerances.
The likelihood uses H1/L1, the analytic acceptance PSD, and a nonzero offset
from the injected draw. Component probes use the first four parameters except
the tail and FFT probes, whose inputs are recorded in the diagnostic file.

| Scalar objective | Enzyme vs ForwardDiff relative norm error | Best Enzyme vs five-point FD error |
| --- | ---: | ---: |
| Calibration dSO | 1.21e-15 | 8.90e-14 |
| Hamiltonian | 5.96e-16 | 1.74e-12 |
| Full multipolar flux | 6.52e-16 | 3.34e-12 |
| Complex tail | 6.04e-16 | 4.72e-11 |

Finite differences use five logarithmically spaced absolute steps from `1e-2`
through `1e-5`, applied separately to each input. The table reports the best
agreement across those steps; all step errors and gradient vectors are saved.
The tail test exercises native complex loggamma successfully.

The following failures produce no Enzyme gradient, so a numerical gradient error
cannot be assigned:

- Initial conditions: `EnzymeNoDerivativeError` for LAPACK `dgesdd_64_`, reached
  through the Newton solver's `cond(J)` check. Runtime activity has the same failure.
- Local FFT adapter: `EnzymeNoDerivativeError` for FFTW
  `fftw_plan_guru64_dft_r2c` during plan creation.
- Scalar continuous timing root: `IllegalTypeAnalysisException` at the bracket
  sign comparison in `roots.jl`. Disabling strict aliasing has the same failure.
- Evolution, complete waveform, and H1/L1 offset likelihood: `IllegalTypeAnalysisException` at
  `sol.retcode == SciMLBase.ReturnCode.Terminated` in `dynamics.jl`.

For comparison, the complete waveform ForwardDiff gradient reaches `1.37e-6`
relative error against finite differences, and the offset likelihood reaches
`2.47e-7`, at their best tested steps. Neither failed Enzyme call returns a
gradient to compare.

These are the first observed blockers, not proof that subsequent operations have
working adjoints. In particular, implicit-root and continuous-time partials are
currently attached through ForwardDiff-specific methods. Correct reverse support
needs implicit pullbacks and validation after compilation barriers are resolved.
No custom reverse rules or runtime model changes were introduced for these probes.

## Reproduction

```sh
julia +1.12 --project=integration/GWSEOBNR -e 'using Pkg; Pkg.instantiate()'
JULIA='julia +1.12' python3 integration/GWSEOBNR/probe_enzyme.py
```

The Python driver needs Python 3.11+ and uses only the standard library. Each
probe runs in a separate Julia process with a 180-second timeout. Failed probes
cause a nonzero driver exit status; they are recorded rather than asserted as
expected successes. Full errors, gradients, finite-difference comparisons, and
logs are written under gitignored `reference/seobnr/enzyme/`. A compact measured
record is committed in `packages/GWSEOBNR/diagnostics/enzyme-2026-10-19.toml`.
Only the standalone integration environment adds Enzyme; package dependencies
and the advertised differentiation capability are unchanged.
