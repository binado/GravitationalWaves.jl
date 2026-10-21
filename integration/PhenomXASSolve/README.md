# IMRPhenomXAS collocation solve prototype

An isolated experiment for the three phase coefficient systems used by IMRPhenomXAS:
4×4 inspiral, 5×5 intermediate, and 5×5 ringdown. This directory preserves the solve
investigation; the complete waveform now lives in `packages/FrequencyDomainWaveforms` and its
ForwardDiff/Mooncake acceptance tests in `integration/FrequencyDomainWaveforms`. Production now
uses StaticArrays and pivoted solves. The Reactant/Enzyme dependencies and explicit
pullback below belong to this historical experiment, outside the supported waveform AD contract.

`collocation.jl` reconstructs the actual solve matrices with StaticArrays and
uses pivoted LU through `A \ b`. It also assembles dense columns directly for
Reactant and solves those systems with `lu(A) \ b`. Inputs are the dimensionless
MECO, ISCO, ringdown and damping frequencies, followed by the assembled collocation RHS.
All of these are active inputs; this experiment does not yet differentiate the
upstream remnant fits or fitted collocation values with respect to masses/spins.

`fixtures.toml` contains 18 systems from six mass/spin configurations, including
equal masses, unequal spins, positive and negative aligned spins, and a 100:1
mass ratio. The reference formulas come from
[Ripple revision 58c5302](https://github.com/tedwards2412/ripple/tree/58c53024718697e0430b666ef5dd40ff3c2e029d/src/ripplegw/waveforms),
evaluated using NumPy float64, not JAX. These fixtures establish solve-level
agreement with those formulas, not complete waveform agreement with LAL.
The unit conversion cancels from this calculation and is set to one during
fixture generation. The intermediate inspiral endpoint derivative is obtained
by complex step solely for preparing its fixture RHS.

From the repository root, with Julia 1.12:

```sh
julia --project=integration/PhenomXASSolve -e 'using Pkg; Pkg.instantiate()'
julia --project=integration/PhenomXASSolve integration/PhenomXASSolve/runtests.jl
julia --project=integration/PhenomXASSolve integration/PhenomXASSolve/test_reverse.jl
julia --project=integration/PhenomXASSolve integration/PhenomXASSolve/test_reactant.jl --values-only
julia --project=integration/PhenomXASSolve integration/PhenomXASSolve/test_reactant_pullback.jl
```

The first script checks matrix reconstruction, captured solutions, normwise
backward errors, a 256-bit BigFloat solve, and full ForwardDiff Jacobians against
finite differences and the implicit equation `dc = A \ (db - dA*c)`. Inputs and
coefficients are scaled for comparison; no change of basis or equilibration is
applied to the actual solve. A second draw perturbs both frequencies and RHS.

The reverse script compares Enzyme and Mooncake gradients to ForwardDiff.
The Reactant value script compiles one executable per region and reuses it
across all six fixture draws and their perturbations. These include different
LU pivot patterns for the intermediate system.

`test_reactant_pullback.jl` evaluates the real-valued implicit solve pullback:

```text
c = A \ b
lambda = transpose(A) \ cbar
Abar = -lambda * transpose(c)
bbar = lambda
```

It then uses Enzyme to differentiate only the assembly, with `Abar` and `bbar`
as fixed cotangents for that invocation. This includes the matrix dependence.
An independent host check uses ForwardDiff on the assembly and compares the
explicit chain rule with ForwardDiff through the complete StaticArrays solve.
It tests an explicit first-order gradient pipeline, not a registered AD rule for
arbitrary calls to `\`. Registering/integrating a rule with Reactant is separate
work; higher-order derivatives are not established. The prototype factors the
transpose separately to avoid dependence on backend transpose-factorization
wrappers. It never forms a matrix inverse.

For diagnostics, omit `--values-only` to attempt automatic Reactant AD through
LU itself. Add `--probe-static` to also attempt the StaticArrays and
static-to-dense conversion paths. Unsupported paths print their errors and
return a nonzero exit code; a caught compilation failure does not count as a
passing check.

## Results and interpretation

Tested with Julia 1.12.0, StaticArrays 1.9.23, ForwardDiff 1.4.6,
FiniteDifferences 0.12.34, Enzyme 0.13.216, Mooncake 0.5.54, and Reactant
0.2.295 on CPU. The ignored local Manifest records the resolved environment.

- 414 matrix/solution/ForwardDiff checks pass, including the dense assembly.
- 144 Enzyme/Mooncake reverse-mode checks pass against ForwardDiff.
- 72 Reactant compiled-value checks pass with direct dense assembly.
- 144 explicit pullback checks pass: the host chain rule and compiled Reactant
  gradients agree with ForwardDiff through the complete solve. The same
  executable per region serves all six configurations and their perturbations.
- Direct automatic Reactant AD of dense LU fails: the installed stack reports
  missing adjoints for `enzymexla.linalg.lu` and (for general backslash's
  structural branches) `stablehlo.triangular_solve`.
- The tested StaticArrays tracing/conversion paths also fail, including
  static-matrix indexing recursion. Static size alone does not guarantee
  Reactant support.

The largest matrix condition number is about `3.37e8`. For the StaticArrays
path, the maximum coefficient error against the 256-bit solve, divided by each
coefficient's reference scale, is `1.72e-12`; the maximum relative Jacobian
error against the implicit equation is `1.03e-12`. Maximum normwise backward
error is `1.84e-17`. These are results on these fixtures, not a guarantee over
the full model's calibration domain. No performance comparison was made.

The experiment established StaticArrays as a viable CPU starting point. Direct dense assembly plus the
implicit solve pullback gives a working first-order compiled gradient pipeline
in this prototype. Production now targets ForwardDiff and Mooncake only, so it uses
ordinary StaticArrays pivoted solves without a custom pullback. Automatic AD through
the LU primitive was unsupported in the Reactant version tested here.

To regenerate the fixtures, download `IMRPhenomXAS.py` and `IMRPhenomX_utils.py`
from the pinned revision into a directory and run:

```sh
python3 integration/PhenomXASSolve/generate_fixtures.py /path/to/source/directory
```

This requires NumPy. The generator verifies the source SHA-256 hashes before
evaluating definitions; Julia tests require no Python or source download.

The prototype follows the repository's GPL-2.0-or-later license for native
waveform code derived from LAL.
