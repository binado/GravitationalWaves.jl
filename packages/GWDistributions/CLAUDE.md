# Repository Guidelines

GWDistributions.jl: population-distribution building blocks for gravitational-wave source
modeling (Distributions.jl-based), a port of `astrogwb.distributions` (Python/NumPyro): interpolated
distributions, redshift distributions (Madau–Dickinson, optionally time-delayed), delay-time,
max-of-two and orientation distributions, the default BBH mass law and aligned spins.

- Layout: `src/base/` (power laws, `InterpolatedDistribution`), `src/mass/` (BBH mass law,
  `MaxOfTwo*`), `src/spins/`, `src/redshift/` (`cosmology.jl`, `base.jl`, `time_delayed.jl`,
  `madau_dickinson.jl`), `src/delay.jl`, `src/orientation.jl`, `src/utils.jl`.
- Cosmology comes from BackgroundCosmology.jl (unregistered, so it is a `[sources]` URL entry
  in `Project.toml`). `params` is a `NamedTuple` whose keys are BackgroundCosmology's
  (`H0`, `Ωm`, optional `w0`, `wa`) plus the model keys, so `cosmology(params)` builds the model
  straight from it. Do not add other cosmology dependencies.
- Naming follows Distributions.jl (no `Distribution` suffix where it would use none);
  `astrogwb` is the behavioral reference. `test/test_astrogwb_reference.jl` pins values
  printed by astrogwb (JAX, x64); regenerate them from astrogwb's `.venv` if semantics change.
- Quadrature is fixed-order Gauss–Legendre (FastGaussQuadrature), never adaptive, and iterative
  solves use a fixed iteration count, so ForwardDiff duals propagate. Adaptive `QuadGK` is a
  test-only reference.
- Format with Runic (`just install-runic` once, then `just fmt`; CI runs `just fmt-check`).
- Test with `just test` (Julia 1.12; includes Aqua, and JET on 1.12 only).
- Conventional Commits.
