# AbstractWaveformsTests agent guidance

AbstractWaveformsTests owns reusable contract, reference-generator, and gradient suites for AbstractWaveforms backends.

- Keep shared contracts in `src/contract.jl`, exported for backend tests to call;
  they must remain light enough for ordinary `Pkg.test()`.
- Keep heavy reference-generator and AD suites in tagged package test items,
  explicitly selected locally and outside default CI; preserve time/phase-maximized
  mismatch semantics based on GWMatchedFiltering and explicit TD window/padding preparation.
- Route every gradient backend and finite-difference call through
  DifferentiationInterface. Compare backends with each other and FiniteDifferences;
  preserve the Mooncake extension's treatment of fixed objective configuration.
- Use `data/parameters.csv` as the shared physical case source (solar masses, Mpc,
  radians, LAL-source-frame Cartesian spins). Select cases by generator capability;
  do not silently project unsupported precession onto aligned-spin models.
- Keep grids, PSD metrics, windows, model data, and other preparation outside
  differentiation; test parameter gradients through the supplied adapters.
- Do not add JET checks here: contract checks intentionally call method-less
  interface functions. Test consuming generators when changing shared contracts.

See the [README](README.md) and [testing guide](../../README.md#development-and-tests)
for suite usage, parameter conventions, and acceptance recipes.
