# NonlinearSolve initial-condition prototype

`nonlinearsolve_initial.jl` is a standalone candidate for the initial-condition
solve. Production waveforms still use the original solver. No NonlinearSolve or
SciMLSensitivity dependency is added to the workspace package.

The prototype solves the same conservative frequency/radial-force equations and
dissipative radial-momentum equation as GWSEOBNR. It combines them into a single
three-variable system `(r, p_phi, p_rstar)`. The conservative equations are
independent of radial momentum, so their roots agree with the original sequential
construction. It uses `NewtonRaphson` with ForwardDiff Jacobians and absolute/
relative tolerances of `1e-14`. The start frequency is fixed at 16 Hz, matching
the production generator's default `0.8 * minimum_frequency`.

The configured implicit adjoint is:

```julia
SteadyStateAdjoint(; autojacvec = false)
```

Enzyme reverses through SciML's native solve rule. That rule uses dense
ForwardDiff residual Jacobians and a transpose linear solve for the implicit
pullback. It does not differentiate the Newton iterations or route the whole
initial-condition function through ForwardDiff. These are first-order, mixed
forward/reverse derivatives.

The solution validator remains active in the primal and inactive for Enzyme. It
checks the return code, residual, Jacobian conditioning, positive radius/angular
momentum, and the same negative radial-momentum branch used by the original
bracket. Failed solves are rejected.

## Reproduction

```sh
JULIA='julia +1.12' just test-integration-seobnr-initial
```

Or use the prototype directly:

```julia
include("integration/GWSEOBNR/nonlinearsolve_initial.jl")
using .NonlinearInitialPrototype
using Enzyme
P = NonlinearInitialPrototype
x = [40.0, 20.0, 0.4, -0.3] # solar masses and signed aligned spins
u = P.initial_conditions(x)  # (r, phi, p_rstar, p_phi)
gradient = zeros(4)
Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(P.initial_loss), Enzyme.Active,
    Enzyme.Duplicated(copy(x), gradient))
```

Use `--project=integration/GWSEOBNR`. The runner checks the nine combinations of
total masses 40, 60, and 100 solar masses with `(q, chi1, chi2)` equal to
`(1,0,0)`, `(2,0.4,-0.3)`, and `(4,-0.6,0.6)`. Every initial-state component must
match the original solver to `rtol=1e-9`, `atol=1e-13`, and residuals must be
below `1e-12`. Each nonconstant state's complete mass/spin gradient is checked
against the original solver's ForwardDiff derivatives and five-point finite
differences at five logarithmically spaced steps. Every gradient row must agree
to `1e-3`, including at two adjacent finite-difference steps. Tests use internal
signed-spin/mass continuations at envelope boundaries; they do not extend the
public waveform's domain or independent-mass derivative contract at equal mass.
The module's default configuration and an exhausted iteration budget are also
checked. Unmet gates produce a nonzero exit status.

## Measured results

On Julia 1.12.7, NonlinearSolve 4.32.0, SciMLSensitivity 7.119.13, and Enzyme
0.13.216, the dense-adjoint candidate passes all nine physical cases. Maximum
relative errors are `5.11e-11` for the initial-state vector (`4.98e-10` for an
individual component) and `6.71e-10` for individual
Enzyme gradient rows against the reference ForwardDiff derivatives. The worst
best-step finite-difference row error is `9.90e-8`; every row also has two
adjacent steps with errors at most `3.85e-7`. Full state/Jacobian/error records are committed under
`packages/GWSEOBNR/diagnostics/nonlinearsolve-initial-2026-10-19.toml`.
The module-default gradient check passes, and an exhausted iteration budget
returns `MaxIters` and is rejected explicitly.

The automatic sensitivity selection (`sensealg=nothing`) was also tried. For this
out-of-place residual it selects Zygote and fails on mutation of coefficient
arrays. That failure is recorded rather than classified as a passing test:

```sh
julia +1.12 --project=integration/GWSEOBNR \
    integration/GWSEOBNR/probe_nonlinearsolve.jl default default \
    reference/seobnr/enzyme/nonlinearsolve-default.toml
```

The explicit dense configuration avoids that backend. This validates an
initial-condition candidate, not complete waveform reverse-mode support. ODE
sensitivities and continuous attachment/peak timing still need separate work.
