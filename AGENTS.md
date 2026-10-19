# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## What this is

GravitationalWaves.jl is a monorepo of composable Julia packages for gravitational-wave modeling.
The root is a Julia workspace, not an umbrella package. Packages retain independent APIs,
versions, dependencies, and licenses.

## Repository layout

Julia packages live in `packages/`:

- `GWTimeFrequencyDomains`: detached `TimeDomain`/`FrequencyDomain` metadata, explicit host checks, crop/pad operations over domain-array pairs, and optional FFTW conversion/plans. Numerical evaluation is domain-first and does not repeat compatibility checks. Domains are fixed configuration in the initial compiled contract.
- `GWNoise`: native one-sided PSD/ASD text readers, tabulated interpolation and finite Gaussian noise segments. Independent of LAL/Python; tests cover physical Fourier normalization. PSD estimation and streaming remain out of scope.
- `GWMatchedFiltering`: joint frequency-domain metrics over ordered channels, real/complex inner products, SNRs and fixed-PSD Gaussian likelihood statistics. GWTimeFrequencyDomains/LinearAlgebra/StaticArrays runtime; independent PSDs use diagonal weights, correlated spectra use cached whitening factors per frequency. Shapes are `(frequency, channel)`; layout measurements live in `benchmarks/GWMatchedFiltering`. Native host matrices, vectors, and their views use allocation-free kernels; traced arrays use broadcast/reduction kernels. Heavy AD and waveform/projection composition live in `integration/GWMatchedFiltering`. Only strain parameters, not metric preparation, are in the initial AD contract. Matched-filter SNR assumes a finite nonzero template; run `check_template` explicitly on the host.
- `AbstractGWWaveforms`: dependency-free interface owning `AbstractWaveformGenerator`, the generic `hphc_fd(f, params, generator)`, the accessors (`approximant`, `reference_frequency`, `minimum_frequency`, `maximum_frequency`), capability queries (`spin_support`, `supports_mode_array`, `is_differentiable`), and the shared helpers `check_frequency_grid`/`frequency_bounds`. `test/contract.jl` defines `test_generator_contract`, which both waveform packages `include` from their tests.
- `LALWaveforms`: direct `ccall` bindings to LALSimulation FD waveforms, matching bilby's parameter names, units, conversions, and floating-point behavior.
- `GWParameterConversion`: dependency-free, AD-friendly mass, spin, and frame conversions ported from bilby.
- `GWWaveforms`: pure-Julia, AD-friendly TaylorF2 and IMRPhenomD, with the same API as `LALWaveforms` (generator types `TaylorF2`, `IMRPhenomD`); tests use LAL as the oracle.
- `GWDetectors`: StaticArrays-based, differentiable interferometer geometry, LAL-convention antenna factors, native GPS-to-GMST, arrival delays, and domain-first broadcasted FD projection over uniform integer-indexed frequency bins. Heavy CPU Reactant/Enzyme/Mooncake acceptance lives in `integration/GWDetectors`.
- `GWDistributions`: Distributions.jl-based population models for masses, spins, orientations, delays, and redshifts. Cosmology comes from the external, unregistered `BackgroundCosmology.jl` dependency via a `[sources]` URL.

The waveform packages depend on `GWParameterConversion` through relative `[sources]` paths.
Shared bilby fixture tooling lives in `python/`; generated fixtures go to `reference/` (gitignored).
Waveform scope remains FD only, uniform `AbstractRange` grids starting at 0, detector-frame masses, no tides.
The waveform backends, conversion, detector, series, noise, and matched-filtering packages use GPL-2.0-or-later; `AbstractGWWaveforms` and `GWDistributions` are MIT.
Paths in the architecture section below are relative to `packages/LALWaveforms/` unless noted.

## Commands

LAL is not bundled; tests need `LALWAVEFORMS_LIBDIR` pointing at the `liblal`/`liblalsimulation` directory (see README).

```sh
# One-time: generate bilby reference fixtures (not committed) and locate LAL from the lalsuite wheel
uv sync --project python
uv run --project python python python/generate.py
uv run --project python python python/generate_detectors.py
export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)

julia --project=packages/LALWaveforms -e 'using Pkg; Pkg.test()'     # full suite (Aqua, JET, reference, etc.)
julia --project=packages/GWDistributions -e 'using Pkg; Pkg.test()'
julia --project=packages/GWParameterConversion -e 'using Pkg; Pkg.test()'
julia --project=packages/GWDetectors -e 'using Pkg; Pkg.test()'
julia --project=packages/GWWaveforms -e 'using Pkg; Pkg.test()'   # needs LALWAVEFORMS_LIBDIR for the LAL comparison
runic --inplace .                              # format (Runic; CI enforces it)
just test-integration                         # heavy AD integration suites, local only
just test-integration-conversion              # Reactant CPU + Enzyme, no LAL/Python needed
just test-integration-waveforms               # ForwardDiff + Enzyme + Mooncake + Reactant CPU (TaylorF2)
just test-integration-detectors               # detector AD + Reactant CPU
just test-integration-inference               # metric AD + waveform/projection/likelihood
just benchmark-inference                      # compare frequency/channel storage layouts
```

Single test file: `julia --project=packages/LALWaveforms -e 'using Test, LALWaveforms; include("packages/LALWaveforms/test/test_fd.jl")'` (the `test/` extras like JSON/NPZ/JET are only resolved under `Pkg.test()`, so for reference/quality tests use `Pkg.test()` or activate `test/`'s target env).

Without fixtures the reference tests warn and skip; with `CI=true` they fail.

## Architecture

`packages/LALWaveforms/src/LALWaveforms.jl` includes files in dependency order: `types` → `library` → `errors` → `generator` → `parameters` → `conversion` → `fd`.

- **`library.jl`** – lazy loading of `liblal`/`liblalsimulation`. `__init__` must never throw (the package must load and pass Aqua without LAL), so libraries are dlopen'd on first `lal()` call. Symbols are listed in `LIBLAL_SYMBOLS`/`LIBSIM_SYMBOLS` and expanded via `@eval` into the `LALFunctions` struct with one `Ptr{Cvoid}` field each, so `lal().XLALFoo` is a type-stable field read. **To call a new XLAL function, add its symbol to these tuples first.** The libdir comes from `LALWAVEFORMS_LIBDIR` (wins) or a saved Preferences entry.
- **`generator.jl`** – `LALWaveformGenerator <: AbstractWaveformGenerator` holds everything constant across draws (approximant, f_min/f_ref/f_max, PN orders) plus an owning `LALDictHandle` (finalizer calls `XLALDestroyDict`). Mirrors bilby's `set_waveform_dictionary`. Its `spin_support` comes from `XLALSimInspiralGetSpinSupportFromApproximant` once at construction. Callers must `GC.@preserve` the generator across `ccall`s that use `dict.ptr`.
- **`parameters.jl` / `conversion.jl`** – `bilby_parameters` validates the NamedTuple (unknown keys throw), derives `mass_1`/`mass_2` from other mass parameters via `GWParameterConversion.convert_parameters`, and returns a plain `Float64` NamedTuple (no parameter struct); `bilby_to_lalsimulation_spins` ports bilby's conversion, short-circuiting to aligned spins and otherwise calling `XLALSimInspiralTransformPrecessingNewInitialConditions`.
- **`fd.jl`** – `hphc_fd` dispatches to `XLALSimInspiralChooseFDWaveform` (native FD approximants) or `XLALSimInspiralFD` (time-domain approximants transformed to FD, followed by bilby's time shift). Output is copied into Julia arrays, zero-masked outside `[f_min, f_max]`, and LAL series are always destroyed in `finally`.
- **`errors.jl`** – XLAL calls use a silent error handler; each call does `clear_errno!` → `ccall` → `check(rc, name, fns)`, which raises `LALError`.

## GWWaveforms

`packages/GWWaveforms/src/` is generic over `Real` (ForwardDiff, Enzyme and Mooncake all differentiate `hphc_fd`; default tests check ForwardDiff against finite differences, and the local integration suite checks all three backends). `taylorf2.jl` ports `XLALSimInspiralPNPhasing_F2` and `XLALSimInspiralTaylorF2`; `generators.jl` defines the generator types (`TaylorF2` and `IMRPhenomD`, each with only the fields it uses); `phenomd.jl` ports `LALSimIMRPhenomD*.c`, with the calibrated fits in `phenomd_fits.jl` and the QNM tables in `phenomd_qnm.jl` (transcribed mechanically from LAL's C, keeping its operation order). LAL wins over ripple or the papers wherever they differ. Keep the kernels free of `try`/`catch` and runtime dispatch (JET and the AD backends need that); validation runs before them and may throw. `hphc_fd` does not check parameter values; `validate_parameters` (host, once) does, so `hphc_fd` has no branches on parameters. `hphc_fd(f, params, generator)` compiles under Reactant CPU for both models via `taylorf2_fd_slice` / `phenomd_fd_slice`, which evaluate a static index slice of `f` (`taylorf2_slice`) and pad with zeros; `taylorf2_fd` and `phenomd_fd` are the host reference loops, kept for the tests. `maximum_frequency` must be positive (the mass-dependent default cutoff cannot be compiled). PhenomD's only parameter-dependent edge, `fCut = 0.2 / M`, is a dynamic mask (`i < floor(min(fCut, f_max) / Δf)`, masked samples are evaluated at an in-band `Mf` and zeroed); the region blending is nested `ifelse` over pieces that are all finite for `Mf > 0`, `spline_eval` selects its segment with a one-hot sum (Reactant rejects scalar indexing with a traced index), and the `eta` and final-spin checks live in `check_waveform_inputs`, not in `phenomd_model`. Avoid `findfirst`/`(i:j) .* x` on float ranges in kernels (Enzyme cannot differentiate `floatrange`). Both models reproduce `LALWaveforms` to about 1e-12 relative (largest sample), bar libm differences on Linux.

## GWDistributions

Follow the package-local guidance in `packages/GWDistributions/AGENTS.md`.
`astrogwb.distributions` is the behavioral reference; existing test values are pinned from JAX x64.
Keep fixed-order Gauss–Legendre quadrature and fixed iteration counts so ForwardDiff duals propagate.
Adaptive QuadGK is test-only. Keep BackgroundCosmology external; do not add other cosmology dependencies.

## Fidelity gotchas

- Match bilby's floating-point operation order (e.g. `((2π) * dt) * f` in the time shift) and compare with `Float64(π)`, not the `Irrational`, when mirroring `np.pi` comparisons.
- Julia loads `libopenlibm` globally, so on Linux LAL's libm calls hit openlibm instead of glibc; results differ at the ulp level from Python (up to ~1e-6 for IMRPhenomT). `python/cases.json` therefore carries a per-case `rtol`. macOS is bit-identical. See README for the `LD_PRELOAD` workaround.
- `python/generate.py` produces the fixtures from bilby; new reference cases go in `cases.json`.

## Git

Conventional Commits; work on a non-`main` branch and open PRs.

`GWTimeFrequencyDomains` changes also test `GWNoise`, `GWDetectors`, and `GWMatchedFiltering`.

PR CI tests affected packages on Julia 1.12, including downstream test dependencies:
`GWParameterConversion` changes also test both waveform packages, `LALWaveforms`
changes also test `GWWaveforms`, and `AbstractGWWaveforms` changes test both waveform packages. Python tooling changes test all three conversion/waveform
packages, `GWDetectors`, and `GWMatchedFiltering`. Detector changes select `GWDetectors`. Inference changes select `GWMatchedFiltering`. CI configuration, root Julia environment files, and unfamiliar paths test all
packages. Only changes confined to root `README.md`, `AGENTS.md`, `CLAUDE.md`, and
`LICENSE` skip package tests; package-local documentation still selects its package.
Formatting runs on every workflow, coverage collection is disabled, and newer PR commits
cancel superseded runs. `CI passed` aggregates detection, test, and formatting results.

Heavy AD checks live in standalone `integration/` environments, outside the root workspace
and default test targets. Run `just test-integration` locally before submitting changes
to numerical kernels or the AD integrations; individual package recipes are listed above.
CI does not run these integration suites or install Enzyme, Reactant, or Mooncake.

Pushes to `main` and manual CI runs test all nine packages on Julia 1.12 and latest stable.
For a full pre-merge run, use `gh workflow run CI.yml --ref <branch>`.
Selection and status checks can be tested locally with
`python3 -m unittest discover -s .github/ci -v`.
