# GWSEOBNR (experimental)

Native Julia aligned-spin SEOBNRv5HM prototype, tied to pySEOBNR **0.3.7**.
**Only (2,±2) waveform output is implemented.** The radiation-reaction flux
retains all upstream multipoles through ℓ=8.

This is an unvalidated research prototype. `is_differentiable` returns **false**
until every primal, derivative, convergence, and composition acceptance gate
passes. ForwardDiff CPU first derivatives are the intended contract; no Hessian,
reverse-mode, Reactant, GPU, precession, eccentricity, or performance guarantee.

Inputs use the shared bilby parameter names and GWParameterConversion mass
conversions. Detector-frame total masses are restricted to 40–100 solar masses,
ordered mass ratios mass_1/mass_2 to 1–4, spin magnitudes to 0–0.6 with aligned
or anti-aligned tilts, positive distances, and zero eccentricity. Fixed nonzero
spin magnitudes can be differentiated with the tilt held at 0 or π.

```
using GWSEOBNR
f = 0:0.0625:2048
p = (mass_1=40., mass_2=20., luminosity_distance=500.,
     theta_jn=0.7, phase=0.3, a_1=0.4, a_2=0.3, tilt_2=π)
g = SEOBNRv5HM()
validate_parameters(f, p, g)
h = hphc_fd(f, p, g)
```

Defaults: minimum/reference frequency 20 Hz (initially required equal),
maximum frequency `nothing` (last grid bin), relative/absolute solve tolerances
1e-11/1e-12, and post-merger duration 1 second. Keep configuration, sampling,
and tolerances fixed during differentiation.

The waveform starts evolving at 0.8 f_min. Orbital phase is referenced to f_ref;
time zero is the continuous amplitude maximum. For K frequency samples,
N=2(K−1), T=1/Δf, and Δt=T/N. Samples cover [−T+t_post,t_post).
A quintic onset spans instantaneous frequency 0.8 f_min to f_min; a final
0.25-second taper starts only after at least 30 QNM damping times. Inadequate
segments, failed roots/solves, missing ISCOs, and attachment clamping requirements
raise explicit errors. The physical spectrum is Δt*rfft with the start-time
phase correction. Shared FFT APIs are unchanged.

The Hamiltonian uses the upstream analytic gradient and Hessian. Initial
conditions use NonlinearSolve's NewtonRaphson on the three coupled conservative/
dissipative equations, with residual, conditioning, and physical-branch checks.
ForwardDiff uses NonlinearSolve's implicit derivatives without requiring
SciMLSensitivity. Loading optional `SciMLSensitivity` selects the tested
`SteadyStateAdjoint(autojacvec=false)`, using dense ForwardDiff residual Jacobians
inside the native Enzyme implicit pullback. Continuous timing partials use
implicit differentiation through ForwardDiff; reverse timing pullbacks are not
implemented. Evolution uses Vern9 dense output. The complex
tail and FFT have package-local ForwardDiff adapters. Timing and matching use
continuous interpolants, with fixed five-point time derivatives for matching.
Physical piecewise formulas are retained.

Loading EnzymeCore enables an optional reverse rule for the local real FFT on
dense Float32/Float64 vectors (one reverse lane). Planning uses scratch storage;
the pullback applies the constant plan's adjoint. This component rule does not
establish reverse-mode waveform support: continuous timing pullbacks remain
outstanding and the complete reverse path remains unvalidated.
Discrete solver-status validation is inactive under Enzyme while retaining its
primal error checks; trajectory values and continuous timing remain active.
Run `just test-integration-seobnr-fft` for direct-DFT and gradient checks.
Run `just test-integration-seobnr-initial` for nine-case production initial-state
and Enzyme-gradient comparisons against the frozen sequential solver.

## Current numerical and reverse-mode errors

Measured on 2026-10-19 with Julia 1.12.7, Enzyme 0.13.216, NonlinearSolve
4.32.0, and SciMLSensitivity 7.119.13. No complete Enzyme waveform gradient has
been returned, so there is no waveform reverse-gradient numerical error to report.

- **Primal/ForwardDiff acceptance:** 130 checks pass and 13 fail. Nine failures
  concern NQC matching coefficients; four concern the unweighted spectral L²
  error limit (`1e-3`). All tested ForwardDiff derivative checks, noise-weighted
  mismatches, convergence checks, and offset likelihood gradients pass. See
  [ACCEPTANCE.md](ACCEPTANCE.md) for the cases and diagnostics. The failed gates
  remain required; `is_differentiable=false`.
- **Production reverse path:** the initial-condition gradients pass all nine
  tested cases, but the complete waveform fails in automatic ODE sensitivity
  selection with `ReverseDiffVJP{false} does not support Functors.jl parameter
  structs`. The ODE parameters are a NamedTuple containing Hamiltonian data and
  precomputed flux arrays. Enzyme is the outer engine; SciMLSensitivity selects
  ReverseDiff internally for this ODE adjoint. See the committed
  [production diagnostics](diagnostics/nonlinearsolve-production-2026-10-19.toml).
- **Experimental SciMLStructures parameters:** an integration-only wrapper
  exposes all 7,297 real coefficient values as tunable, including the real and
  imaginary parts of the multipole prefixes. This gets the full waveform probe
  past the parameter-container error. Its next observed failure is
  `IllegalTypeAnalysisException` in `bracket_root`, at
  `signbit(fa) != signbit(fb)` (`src/roots.jl`). Correct reverse timing requires
  implicit pullbacks; marking the root inactive would discard physical
  derivatives. The wrapper is not installed in production.
- **Packing experiment:** the first short-ODE attempt with the wrapper failed
  with `Need an adjoint for constructor Base.Generator` during SciMLSensitivity's
  Zygote reconstruction of structured parameter cotangents. The revised wrapper
  uses explicit concatenation and gets past that error, but the short ODE now
  fails with `FieldError: type Active has no field dval` in
  `DiffEqBaseEnzymeExt.reverse` (`solve_up`): the initial `SVector` is annotated
  `Active`, while the rule accesses `.dval` on it. No short-ODE reverse gradient is
  returned. Direct RHS Enzyme gradients agree with ForwardDiff to `1.9e-15`;
  this component result does not validate the ODE adjoint or complete waveform.

See [the SciMLStructures experiment](../../integration/GWSEOBNR/SCIMLSTRUCTURES.md),
[committed experiment diagnostics](diagnostics/scimlstructures-2026-10-19.toml),
and [the Enzyme investigation](../../integration/GWSEOBNR/ENZYME.md) for
reproduction commands and the distinction between previous and current blockers.
The initial-condition SVD and FFT-planning errors from earlier probes have been
resolved for their tested components; they are not the current production blocker.

Small committed component/QNM snapshots run without Python. The separate locked
oracle environment is documented in [python/seobnr](../../python/seobnr/README.md).
Full waveform/gradient acceptance belongs to `integration/GWSEOBNR`; large
fixtures are gitignored. See `ACCEPTANCE.md` for measured status and limitations.

GPL-3.0-or-later; see [NOTICE.md](NOTICE.md), [LICENSE](LICENSE), and
[provenance.toml](provenance.toml). Existing workspace package licenses are unchanged.
