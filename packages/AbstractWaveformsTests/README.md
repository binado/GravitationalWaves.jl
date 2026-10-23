# AbstractWaveformsTests

Shared test suites for any generator implementing the
[AbstractWaveforms](../AbstractWaveforms) interface. Backends keep their test suites in their own
packages and depend on this one. MIT licensed.

## Suites

- **Contract**: `test_generator_contract`, `test_amplitude_phase_fd_contract`,
  `test_modes_td_contract`, `test_modes_fd_contract`. Light enough for `Pkg.test()`.
- **Reference**: `test_reference(target, reference, cases; grid, max_mismatch, ...)` compares a
  target generator with a trusted reference generator (`LALWaveforms`, the PyCall `gwsurrogate`
  adapter, ...) by `1 - overlap`. `mismatch(a, b, metric; maximize)` maximizes over a common time
  and phase shift (`:time_phase`, the default), a phase (`:phase`) or nothing (`:none`). Inner
  products come from GWMatchedFiltering; `flat_metric` gives a unit PSD. TD generators go through
  a caller-supplied `prepare` (window, padding) and the physical `rfft`, so FD and TD generators
  share one comparison. Load an FFT backend (for example FFTW) for TD references.
- **Gradients**: `test_gradients`, `test_snr_gradients`, `test_phase_gradients` and
  `test_waveform_jacobian` need no reference generator. They compare the given AD backends with each
  other and with `FiniteDifferences`, always through DifferentiationInterface. Load Mooncake to enable
  the extension that marks `FixedObjective` non-differentiable.

## Cases

`data/parameters.csv` holds physical parameters (solar masses, Mpc, radians, Cartesian spins in the LAL
source frame). `applicable_cases(generators...)` keeps the rows every generator supports according to
`spin_support`; `fd_parameters` and `surrogate_parameters` convert a row to a generator's parameters,
and `scaled_coordinates` with `fd_adapter`/`surrogate_adapter` give dimensionless AD coordinates.

## Scope

The reference and AD suites are heavy and meant to be run locally (see the `just test-reference-*` and
`just test-integration-*` recipes in the repository root), not in CI. There is no JET test here: the
contract checks call interface functions that have no methods of their own.
