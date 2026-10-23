# SurrogatePyCallWaveforms agent guidance

SurrogatePyCallWaveforms adapts upstream gwsurrogate through PyCall as the reference generator for native surrogates.

- Implement the shared AbstractWaveforms TD mode contract; leave shared harmonic
  assembly to AbstractWaveforms and preserve physical strain units, ordered
  signed modes, and viewing convention `(inclination, π/2 - phi_ref)`.
- Keep this package outside the root workspace and CI matrix: its tests run
  locally and require PyCall built against a Python containing gwsurrogate.
- Use the pinned `python/surrogates` environment for reference acceptance. Set
  `PYTHON` before instantiating/building PyCall, and rebuild it after changing
  interpreters; changing the variable at runtime does not switch an existing build.
- Reuse model HDF5 files from the native Scratch cache or a supplied data directory;
  do not download a second copy. Python evaluation is not differentiable.
- Run the package's own `Pkg.test()` after Python setup for contract and validation
  checks; run `just test-reference-surrogates` for native-vs-upstream acceptance.
  The root recipe configures Python for the standalone reference environment.

See the [README](README.md) for adapter setup and
[reference-suite guide](../../integration/SurrogateWaveforms/README.md) for the current pinned environment.
