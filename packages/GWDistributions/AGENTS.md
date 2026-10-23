# GWDistributions agent guidance

GWDistributions provides Distributions.jl population models following astrogwb.distributions behavior.

- Treat `astrogwb.distributions` as the behavioral reference.
  `test/test_astrogwb_reference.jl` pins JAX x64 values; regenerate from astrogwb's
  environment when semantics change.
- Follow Distributions.jl naming, omitting a `Distribution` suffix where it does.
- Keep BackgroundCosmology external via the `[sources]` URL in `Project.toml`;
  do not add other cosmology dependencies.
- Preserve `params` as a NamedTuple combining BackgroundCosmology keys (`H0`, `Ωm`,
  optional `w0`, `wa`) with model keys, so `cosmology(params)` constructs the model.
- Use fixed-order Gauss–Legendre quadrature and fixed iteration counts so
  ForwardDiff duals propagate. Adaptive QuadGK is a test-only reference.

See the [README](README.md) for distribution conventions and cosmology usage.
