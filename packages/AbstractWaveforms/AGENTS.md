# AbstractWaveforms agent guidance

AbstractWaveforms owns the shared frequency- and time-domain waveform-generator interfaces and harmonic assembly.

- Define shared generics, capability queries, validation/grid helpers, and default
  assembly here; keep backend model physics and restrictions in their packages.
  Frequency-domain support is backend-specific; do not impose FD-only scope on TD generators.
- Coordinate-first polarization calls return ordinary `(sample, 2)` plus/cross
  matrices. Series prototypes return series with `(:plus, :cross)` channels:
  FD output uses epoch zero; TD output preserves peak-relative epoch and spacing.
- TD `dt` and uniform `times` keyword conveniences return TimeSeries. Lower-level
  mode APIs may retain native irregular grids; do not automatically FFT between domains.
- Preserve ordered, unique signed mode labels and `(sample, mode)` storage.
  FD modes provide both positive and negative frequency branches; assembly must
  not infer symmetry between signed modes. Modes exclude observer inclination.
- Keep grids/generators fixed during numerical work and host validation separate
  from backend kernels. Preserve promoted numeric types and finite pole derivatives
  in harmonic assembly.
- Shared contracts live in `AbstractWaveformsTests/src/contract.jl`; use that
  package's exported checks in backend tests, rather than including a local
  `test/contract.jl`. Test downstream generators when changing this interface.

See the [README](README.md) for interface details and
[AbstractWaveformsTests](../AbstractWaveformsTests/README.md) for current contract ownership.
