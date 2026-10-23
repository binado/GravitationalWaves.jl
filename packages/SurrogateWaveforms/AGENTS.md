# SurrogateWaveforms agent guidance

SurrogateWaveforms is a differentiable time-domain port of the NRSur7dq4 family with upstream gwsurrogate as reference.

- Preserve upstream node bookkeeping and frame/mode conventions, but retain the
  corrected natural spline in the final interval of short knot windows; do not
  copy gwsurrogate's isolated-point interval-selection defect.
- Store model data in flat integer-indexed vectors with plain loops, without
  dictionaries, closures, or generators over loaded arrays. Keep test oracles
  SphericalFunctions and DataInterpolations out of runtime dependencies.
- Avoid branches on differentiated values inside recurrences; use `element_type`
  for mixed dual/float spin tuples. Preserve guarded zero-spin norms and use fixed
  `t_ref` for gradients because `t_ref_from_omega` is only C⁰.
- Keep the shared TD contract: signed inertial modes exclude observer angles,
  viewing angles are `(inclination, π/2 - phi_ref)`. Coordinate-first polarizations
  return matrices; series prototypes and uniform `dt`/`times` conveniences return
  TimeSeries. Lower-level modes may use native irregular grids.
- Tests need model data and fixtures. Data download to Scratch on first use or
  come from `GWSURROGATES_DATA_DIR`. From the root, sync `python/surrogates` and run
  `uv run --project python/surrogates python python/generate_surrogates.py --data-dir DIR --download`
  for fixtures under `reference/surrogates/`; missing fixtures fail under `CI=true`.
- Run `just test-integration-surrogates` for ForwardDiff/Mooncake and
  `just test-reference-surrogates` for the configured PyCall upstream comparison.

See the [README](README.md) for AD rules, data attribution, and the intentional spline correction;
see the [integration guide](../../integration/SurrogateWaveforms/README.md) for reference setup.
