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

The initial probes below predate the local FFT rule added in the follow-up
described below. These failures produced no Enzyme gradient, so a numerical
gradient error cannot be assigned:

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
No custom reverse rules or runtime model changes were present in the initial probes.

## FFT reverse-rule follow-up

`GWSEOBNREnzymeCoreExt` now provides a reverse rule for `GWSEOBNR.real_fft` on
dense Float32/Float64 vectors, with one reverse lane. EnzymeCore is an optional
weak dependency; Enzyme remains outside the runtime dependencies. The rule
creates an ESTIMATE plan on scratch storage in the augmented primal, retains it
as constant tape data, and applies its existing AbstractFFTs adjoint to the output
cotangent. Input cotangents accumulate and the consumed output cotangent is
cleared. No methods on FFTW or global inactive annotations are introduced.

This follows the constant-plan/adjoint construction discussed in
[NonlinearSolve #597](https://github.com/SciML/NonlinearSolve.jl/issues/597),
restricted to our out-of-place real-vector wrapper. Planning is not differentiated.
The local rule does not provide whole-waveform reverse-mode support.

```sh
just test-integration-seobnr-fft
```

`test_enzyme_fft.jl` checks lengths 1, 2, 7, 8, 15, and 16 for both scalar types,
with analytic direct-DFT gradients, ForwardDiff, and five-point finite
differences. Linear complex probes cover DC and the highest bin (Nyquist for
even lengths). Repeated FFT calls, existing input cotangents, and the active
physical-spectrum normalization/time-origin correction are also checked.
For the quadratic direct-DFT objective the finite-difference step is `1e-2`;
its five-point formula has no truncation error in exact arithmetic. The physical
time-shift test uses `1e-6` and checks discretization independently.

All 131 assertions pass on Julia 1.12.7 / Enzyme 0.13.216. Maximum relative
analytic direct-DFT gradient errors are `1.27e-15` (Float64) and `7.40e-8`
(Float32). Maximum five-point finite-difference errors are `6.50e-14` and
`7.40e-8`, respectively. The original `sum(abs2, real_fft(x))` probe now passes
with `1.56e-16` relative error against ForwardDiff and best finite-difference
error `8.89e-14`. The committed follow-up measurements are in
`packages/GWSEOBNR/diagnostics/enzyme-fft-2026-10-19.toml`; the initial record is
retained as the pre-rule baseline.

At that stage the complete waveform probe was rerun with the FFT extension loaded. It still
returns `IllegalTypeAnalysisException` at the solver return-code equality check
in `dynamics.jl:91`, before producing a gradient. The local FFT fix therefore
does not remove the independent dynamics/root blockers.

## Solver-status validation follow-up

The two return-code checks now live in `check_evolution_retcode`, which returns
`nothing` and is marked inactive only in the optional EnzymeCore extension.
The check still runs in the primal: failed return codes raise the original
failure error, and a successful solve that reaches its time limit without a
physical termination guard raises the original coverage error. Neither the ODE
solve nor any trajectory values or continuous event times are marked inactive.

```sh
just test-integration-seobnr-validation
```

Three Enzyme assertions check the valid-case gradient and both error branches.
Normal package tests also cover these statuses. Evolution and complete waveform
reverse probes now pass the return-code comparison and fail at LAPACK
`dgesdd_64_`, used by `cond(J)` in the initial-condition Newton solver. No reverse
gradient is returned. The condition check and root mathematics are unchanged;
implicit-root/continuous-time pullbacks remain necessary for correct reverse
support. Measurements are recorded in
`packages/GWSEOBNR/diagnostics/enzyme-validation-2026-10-19.toml`.

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
Only the standalone integration environment adds Enzyme. The package adds the
optional EnzymeCore extension; the advertised differentiation capability is unchanged.

## Production NonlinearSolve initial conditions

Initial conditions now use the tested coupled NonlinearSolve system. Loading
optional SciMLSensitivity selects the dense implicit adjoint, while primal and
ForwardDiff evaluation need no SciMLSensitivity import. The old custom Newton
solver is frozen in test/reference tooling and is absent from the runtime path.

The production primitive passes all nine initial-state/Enzyme-gradient cases,
including per-component finite-difference checks. See [NONLINEARSOLVE.md](NONLINEARSOLVE.md)
and `packages/GWSEOBNR/diagnostics/nonlinearsolve-production-2026-10-19.toml`.

The complete Enzyme waveform probe, retried in a fresh process, now returns an
ODE sensitivity error: automatic selection chooses `ReverseDiffVJP`, which does
not support the Functors parameter structure. No waveform gradient is returned.
This replaces the initial-condition SVD blocker; continuous timing still lacks
reverse pullbacks. Full waveform reverse mode remains unvalidated.

## SciMLStructures ODE parameters

An integration-only `EOBParameters` prototype exposes every precomputed physical
coefficient as a flat real tunable vector. See [SCIMLSTRUCTURES.md](SCIMLSTRUCTURES.md)
for the interface, reproduction commands, and measured limitations. The runtime
parameter container and advertised differentiation support are unchanged.
