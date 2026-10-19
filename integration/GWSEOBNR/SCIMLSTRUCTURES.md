# SciMLStructures ODE parameter experiment

This is an integration-only experiment. Runtime `evolve` still uses its existing
NamedTuple parameters; the package has no new dependency or capability claim.
The implementation follows the
[SciMLStructures interface example](https://docs.sciml.ai/SciMLStructures/stable/example/).

`scimlstructures_parameters.jl` wraps the Hamiltonian scalar coefficients and
all nine flux coefficient arrays in `EOBParameters`. Its `Tunable` portion is a
stable flat real vector: seven Hamiltonian scalars, eight column-major 9×9×11
real arrays, and the real and imaginary parts of the 9×9 multipole prefixes.
All 7,297 entries remain differentiable, including draw-dependent coefficients
precomputed before the solve. No coefficient is hidden in `Constants` or
`Caches`. Calling these arrays caches would be incorrect: the RHS reads their
precomputed values rather than overwriting them on every call.

The immutable replacement contract rebuilds the complete structure with the
new scalar type. Canonicalization allocates a separate vector (`aliases=false`);
replacement does not mutate the original object. Absent portions report false
and canonicalize to `(nothing, nothing, nothing)`. No third-party types or
functions are marked inactive.

Run each experiment in a fresh process:

```sh
julia +1.12 --project=integration/GWSEOBNR -e 'using Pkg; Pkg.instantiate()'
julia +1.12 --project=integration/GWSEOBNR integration/GWSEOBNR/probe_scimlstructures.jl short reference/seobnr/enzyme/structures-short-explicit.toml
julia +1.12 --project=integration/GWSEOBNR integration/GWSEOBNR/probe_scimlstructures.jl waveform reference/seobnr/enzyme/structures-waveform-explicit.toml
```

The short experiment checks 102 assertions, including all nine physical cases:
exact coefficient round trips, replacement ordering, alias behavior, exact RHS
values, and ForwardDiff RHS gradients. A fixed 20-geometric-time-unit evolution
also checks exact primal agreement and ForwardDiff agreement with the original
parameter container. It then compares Enzyme gradients of the RHS and short ODE
objective with ForwardDiff and five-point finite differences at five steps.
Failed reverse probes return a nonzero exit code and retain their diagnostics.
The short ODE includes the production NonlinearSolve initial conditions and
uses the existing automatic ODE sensitivity selection, without an explicit VJP
change. It does not exercise merger timing or termination callbacks.

The waveform experiment reads production `dynamics.jl` and evaluates its
`evolve` definition in the dedicated process with exactly one substitution:
`pars = EOBParameters(c, coeffs)`. It changes no equations, tolerances, callback,
termination conditions, or timing roots. It then runs the existing waveform
probe, including ForwardDiff and finite differences. The `evolution` mode is
also available. The process-local definition is an exploratory harness and is
not loaded by ordinary package tests or waveform evaluation.

An initial canonicalization used a splatted generator expression. The short
ODE reached SciMLSensitivity's Zygote cotangent reconstruction but failed with
`Need an adjoint for constructor Base.Generator`. Explicit concatenation removes
that packing obstacle. The initial raw diagnostics are retained under
`reference/seobnr/enzyme/structures-short.toml`; the revised run uses the
`structures-short-explicit.toml` destination. The first full waveform run got
past parameter-container validation and failed compiling the sign comparison
in `bracket_root`, before returning a gradient. This does not establish that
all later waveform adjoints work.

The revised short ODE fails in `DiffEqBaseEnzymeExt.reverse` for `solve_up` with
`FieldError: type Active has no field dval`. The stack records
`u0::Active{SVector{4,Float64}}`, while the rule accesses `u0.dval` at line 140.
The packing change therefore resolves its initial generator failure but does
not yield an ODE gradient. No alternative initial-state representation or rule
patch was tested in this experiment. The direct RHS Enzyme gradient passes:
relative error `1.90e-15` against ForwardDiff, best five-point error `7.11e-12`.
All 102 structural/primal/ForwardDiff assertions pass. These results support
the container's forward behavior, not complete reverse-mode support.

The revised full waveform confirms the same `bracket_root` sign-comparison
compiler failure. Its scalar primal objective and entire ForwardDiff gradient
are bit-identical to the production-container probe at this draw. No Enzyme
waveform gradient is returned. Compact results are committed in
[`scimlstructures-2026-10-19.toml`](../../packages/GWSEOBNR/diagnostics/scimlstructures-2026-10-19.toml);
full exception stacks remain under gitignored `reference/`. The container
experiment has not been promoted to production.
