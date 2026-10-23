# GWParameterConversion.jl

Dependency-free conversions between gravitational-wave source parameters, following
[bilby](https://bilby-dev.github.io/bilby/)'s names and conventions (`mass_1 ≥ mass_2`,
`mass_ratio = mass_2 / mass_1 ≤ 1`). Numerical data arguments accept `Number` to support traced real scalars,
ForwardDiff, and other AD tools. Physical inputs remain real-valued; complex-valued
physics is not supported.

```julia
using GWParameterConversion

# Key-driven: derive what you need from whatever you have.
convert_parameters((chirp_mass = 28.1, mass_ratio = 0.8, theta_jn = 0.4), (:mass_1, :mass_2))
# -> (chirp_mass = 28.1, mass_ratio = 0.8, theta_jn = 0.4, mass_1 = ..., mass_2 = ...)

mc = chirp_mass(36.0, 29.0)
q = mass_ratio(36.0, 29.0)
m1, m2 = chirp_mass_and_mass_ratio_to_component_masses(mc, q)
```

`convert_parameters(nt, targets)` resolves the route at compile time (pass `targets` as a literal),
never overwrites given keys, and throws `InconsistentParametersError` when redundant inputs disagree.
It covers all mass parameters (and their `_source` versions given a `redshift`), plus `spin_1z`,
`spin_2z`, `chi_eff` and `chi_p` from magnitudes and tilts.

Also provided: `effective_spin`, `effective_precessing_spin`, and `spins_to_lal_frame`, a pure-Julia
equivalent of `XLALSimInspiralTransformPrecessingNewInitialConditions` that agrees with LAL to about
1e-14 (LALWaveforms keeps LAL itself as the test oracle).

The mass-ratio solves from `mass_1`/`mass_2` and the chirp mass use bilby's closed-form cubic root,
which has no real root once `(chirp_mass / mass)^5 > 27/4` (possible for `mass_2`); there bilby
returns `NaN` and this package throws `DomainError`.

## Inference coordinate changes

NamedTuple overloads change between independent mass coordinate pairs:

```julia
x = (; chirp_mass = 28.0, mass_ratio = 0.8)
f = chirp_mass_and_mass_ratio_to_component_masses
masses = f(x) # (; mass_1, mass_2)
component_masses_to_chirp_mass_and_mass_ratio(masses) # recovers x

total = component_masses_to_total_mass_and_mass_ratio(masses)
total_mass_and_mass_ratio_to_component_masses(total) # recovers masses
```

Each NamedTuple must contain exactly the two corresponding coordinate keys, in either
order; outputs use the order shown above. Masses must be finite and positive, with
`mass_1 ≥ mass_2` and `0 < mass_ratio ≤ 1`. Incorrect keys or nonnumeric values throw
`ArgumentError`. Physical domains are caller preconditions and are not checked. Invalid
physical inputs can produce arithmetic errors or invalid results that differ by backend.
Equal masses are supported.
The existing positional methods retain their behavior and return ordinary tuples;
the two new `component_masses_to_*` functions also support that convention.

For analytic Jacobians and inverse registration, install and load the optional
ChangesOfVariables and InverseFunctions packages:

```julia
using ChangesOfVariables: with_logabsdet_jacobian
using InverseFunctions: inverse

y, ladj = with_logabsdet_jacobian(f, x)
inverse(f)(y) # recovers x
_, inverse_ladj = with_logabsdet_jacobian(inverse(f), y) # approximately -ladj
```

`ladj` is the log absolute determinant in the direction of the function called:
`log p_y(y) = log p_x(x) - ladj`. The Jacobian methods cover all four NamedTuple
transformations and are also accessible through Bijectors' re-export of
`with_logabsdet_jacobian`. Each integration loads independently; neither adds a required
dependency to GWParameterConversion. Unlike these coordinate changes,
`convert_parameters` retains inputs and adds derived fields, so it has no Jacobian method.
The transformations and log determinants preserve ForwardDiff duals and BigFloat precision,
including derivatives and nested derivatives at equal masses.

## Skipping the consistency check

`check_consistency=true` is the host default, including `InconsistentParametersError`. With
`false`, redundant-input derivations and comparisons are omitted entirely. Redundant inputs must
therefore be consistent; supplied keys retain precedence and unrelated fields pass through. Keep
target keys, the consistency flag, and tolerance as fixed configuration of differentiated calls.
Enzyme and Reactant compilation were tested earlier and are deferred; see the repository history.

Masses must be positive and finite, with `mass_1 ≥ mass_2`, `0 < mass_ratio ≤ 1`,
and `0 < symmetric_mass_ratio ≤ 1/4`; redshifts must satisfy `1 + redshift > 0`.
Spin magnitudes and angles must describe physical real spins. Frame transforms require
positive reference frequency and nonzero total angular momentum.
The closed-form cubic solves also have the real-root restriction described above.

Derivatives are supported at smooth physical inputs. Do not require derivatives at
nondifferentiable alignment switches, ties in the maximum defining `chi_p`, zero in-plane
spin norms, the equal-mass square-root inverse of symmetric mass ratio, or undefined
frame azimuths. Boundary values are supported separately: aligned results have exactly zero
in-plane components, and line-of-sight degeneracy uses the fixed fallback direction.
The mass-coordinate transforms and their log determinants are smooth at equal masses.

## Tests

Reference tests compare with bilby output generated by `python/generate.py` at the repo root
(see the top-level README); without fixtures they warn and skip, and fail under `CI=true`.

```sh
julia --project=packages/GWParameterConversion -e 'using Pkg; Pkg.test()'
```

The suite includes ForwardDiff checks.
