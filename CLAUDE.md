# Repository Guidelines

GWDistributions.jl: population-distribution building blocks for gravitational-wave source
modeling (Distributions.jl-based): interpolated 1-D distributions, the default BBH mass
law, aligned spins, and the Madau–Dickinson redshift model.

- Layout: `src/base/` (power laws, `Interpolated1DDistribution`), `src/mass/`, `src/spins/`,
  `src/redshift/`, `src/utils.jl`.
- Cosmology-independent by design: redshift distributions take solid-angle-integrated
  `dV_c/dz` **arrays**. Do not depend on BackgroundCosmology anywhere, not even in the test
  target; tests compute a flat-ΛCDM grid inline.
- Quadrature is fixed-order Gauss–Legendre (FastGaussQuadrature), never adaptive, so
  ForwardDiff duals propagate. Adaptive `QuadGK` is a test-only reference.
- Format with Runic (`just install-runic` once, then `just fmt`; CI runs `just fmt-check`).
- Test with `just test` (Julia 1.12; includes Aqua, and JET on 1.12 only).
- Conventional Commits.
