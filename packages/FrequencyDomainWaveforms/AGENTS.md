# FrequencyDomainWaveforms agent guidance

FrequencyDomainWaveforms implements native differentiable FD models with LALSimulation as the numerical reference.

- Prefer LAL over papers or Ripple when they disagree; preserve LAL operation order.
  Scope is uniform `AbstractRange` grids, detector-frame masses, aligned spins,
  quasi-circular binaries, and no tides. Evaluate supplied frequencies directly,
  including nonzero starts; TaylorF2 retains LAL's grid-step-dependent `-1/Δf` shift.
- Run `validate_parameters` on the host; keep grids and generators fixed and
  `maximum_frequency` positive. Kernels remain generic over `Real`, type stable,
  and free of `try`/`catch`; ordinary indexing and parameter branches are allowed.
- Each model owns `amplitude_phase_fd`: prepare coefficients once, then loop over
  one scalar amplitude/phase kernel. Use lazy piecewise branches; zero amplitude
  and phase outside support and preserve the shared polarization convention.
- Preserve PhenomD's exclusive edge `f < Δf * floor(min(0.2 / M, f_max) / Δf)`
  (`M` in seconds), and PhenomXAS's inclusive cutoff at `0.3 / M_sec`.
- Preserve direct component-mass differences, regular near-zero-spin ISCO evaluation,
  and normalized intermediate-amplitude coordinates. Use pivoted StaticArrays
  collocation solves and differentiate both matrix and RHS.
- `Pkg.test()` needs no LAL libraries. Run `just test-reference-waveforms` with
  configured LAL libraries for reference acceptance, and
  `just test-integration-waveforms` for ForwardDiff/Mooncake acceptance.

See the [README](README.md) for model details; use the
[integration guide](../../integration/README.md) for current reference-suite setup.
