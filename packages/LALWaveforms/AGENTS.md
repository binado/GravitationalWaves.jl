# LALWaveforms agent guidance

LALWaveforms binds LALSimulation frequency-domain waveforms while matching bilby's parameters, units, and arithmetic.

- Preserve bilby's operation order, including `((2π) * dt) * f` time shifts and
  `Float64(π)` comparisons. Scope is uniform zero-start FD grids, detector-frame
  masses, and no tides; copy output and zero-mask outside configured bounds.
- Keep library loading lazy: `__init__` must never throw without LAL installed.
  `LALWAVEFORMS_LIBDIR` overrides saved Preferences; load on the first `lal()` call.
- Register new XLAL symbols in `LIBLAL_SYMBOLS`/`LIBSIM_SYMBOLS` in `src/library.jl`
  before calling them, preserving type-stable `LALFunctions` pointer fields.
- Reuse the generator's owning dictionary and `GC.@preserve` the generator around
  every `ccall` using `dict.ptr`. Always destroy returned LAL series in `finally`.
- Preserve silent XLAL error handling: `clear_errno!` → `ccall` →
  `check(rc, name, fns)` raising `LALError`. Keep unknown-parameter validation and
  bilby mass/spin conversion behavior intact.
- For full tests, `uv sync --project python`, then run
  `uv run --project python python python/generate.py` from the repository root.
  Set `LALWAVEFORMS_LIBDIR` to the directory reported by
  `uv run --project python python python/libdir.py`. Add reference cases in
  `python/cases.json`; preserve per-case tolerances for Linux libm differences.
  Missing fixtures warn locally and fail under `CI=true`.

See the [repository README](../../README.md) for library setup and Linux fidelity details.
