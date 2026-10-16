# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## What this is

A Julia package that calls LALSimulation frequency-domain waveforms directly via `ccall` (no PyCall/LAL wrappers), reproducing bilby's `lal_binary_black_hole` parameter names, units and parameter-to-LAL conversions. The goal is bit-level agreement with bilby/Python. Scope is deliberately narrow: FD only, uniform `AbstractRange` grids starting at 0, detector-frame masses, no tidal parameters.

## Commands

LAL is not bundled; tests need `LALWAVEFORMS_LIBDIR` pointing at the `liblal`/`liblalsimulation` directory (see README).

```sh
# One-time: generate bilby reference fixtures (not committed) and locate LAL from the lalsuite wheel
uv sync --project test/python
uv run --project test/python python test/python/generate.py
export LALWAVEFORMS_LIBDIR=$(uv run --project test/python python test/python/libdir.py)

julia --project -e 'using Pkg; Pkg.test()'     # full suite (Aqua, JET, reference, etc.)
runic --inplace .                              # format (Runic; CI enforces it)
```

Single test file: `julia --project -e 'using Test, LALWaveforms; include("test/test_fd.jl")'` (the `test/` extras like JSON/NPZ/JET are only resolved under `Pkg.test()`, so for reference/quality tests use `Pkg.test()` or activate `test/`'s target env).

Without fixtures the reference tests warn and skip; with `CI=true` they fail.

## Architecture

`src/LALWaveforms.jl` includes files in dependency order: `types` → `library` → `errors` → `settings` → `parameters` → `conversion` → `fd`.

- **`library.jl`** – lazy loading of `liblal`/`liblalsimulation`. `__init__` must never throw (the package must load and pass Aqua without LAL), so libraries are dlopen'd on first `lal()` call. Symbols are listed in `LIBLAL_SYMBOLS`/`LIBSIM_SYMBOLS` and expanded via `@eval` into the `LALFunctions` struct with one `Ptr{Cvoid}` field each, so `lal().XLALFoo` is a type-stable field read. **To call a new XLAL function, add its symbol to these tuples first.** The libdir comes from `LALWAVEFORMS_LIBDIR` (wins) or a saved Preferences entry.
- **`settings.jl`** – `WaveformSettings` holds everything constant across draws (approximant, f_min/f_ref/f_max, PN orders) plus an owning `LALDictHandle` (finalizer calls `XLALDestroyDict`). Mirrors bilby's `set_waveform_dictionary`. Callers must `GC.@preserve` the settings across `ccall`s that use `dict.ptr`.
- **`parameters.jl` / `conversion.jl`** – `BBHParameters` validates the NamedTuple (unknown keys throw); `bilby_to_lalsimulation_spins` ports bilby's conversion, short-circuiting to aligned spins and otherwise calling `XLALSimInspiralTransformPrecessingNewInitialConditions`.
- **`fd.jl`** – `hphc_fd` dispatches to `XLALSimInspiralChooseFDWaveform` (native FD approximants) or `XLALSimInspiralFD` (time-domain approximants transformed to FD, followed by bilby's time shift). Output is copied into Julia arrays, zero-masked outside `[f_min, f_max]`, and LAL series are always destroyed in `finally`.
- **`errors.jl`** – XLAL calls use a silent error handler; each call does `clear_errno!` → `ccall` → `check(rc, name, fns)`, which raises `LALError`.

## Fidelity gotchas

- Match bilby's floating-point operation order (e.g. `((2π) * dt) * f` in the time shift) and compare with `Float64(π)`, not the `Irrational`, when mirroring `np.pi` comparisons.
- Julia loads `libopenlibm` globally, so on Linux LAL's libm calls hit openlibm instead of glibc; results differ at the ulp level from Python (up to ~1e-6 for IMRPhenomT). `test/python/cases.json` therefore carries a per-case `rtol`. macOS is bit-identical. See README for the `LD_PRELOAD` workaround.
- `test/python/generate.py` produces the fixtures from bilby; new reference cases go in `cases.json`.

## Git

Conventional Commits; work on a non-`main` branch and open PRs (CI: tests on Julia 1.12 and latest, plus Runic format check).
