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
- `GWWaveforms`: pure-Julia, AD-friendly TaylorF2, IMRPhenomD and IMRPhenomXAS, with the same API as `LALWaveforms` (generator types `TaylorF2`, `IMRPhenomD`, `IMRPhenomXAS`); tests use LAL as the oracle.
- `GWDetectors`: StaticArrays-based, differentiable interferometer geometry, LAL-convention antenna factors, native GPS-to-GMST, arrival delays, and domain-first broadcasted FD projection over uniform integer-indexed frequency bins. Heavy CPU Reactant/Enzyme/Mooncake acceptance lives in `integration/GWDetectors`.
- `GWSurrogates`: MIT, differentiable port of the NRSur7dq4 family from `gwsurrogate` (time-domain). Currently kernels (natural cubic spline, closed-form Wigner-D l <= 5, guarded norms; `SphericalFunctions.jl`/`DataInterpolations.jl` are test oracles) the coorbital waveform surrogate (`load_coorbital`, `coorbital_modes`) and the dynamics surrogate (`load_dynamics`, `dynamics`: RK4 start-up + variable-step AB4 on an `SVector{11}` state, mirroring `gwsurrogate`'s node bookkeeping; `gwsurrogate`'s C++ spline mis-evaluates isolated points in the last interval of a 4-knot window, which this port does not copy); the full pipeline (`load_surrogate`, `inertial_modes`: dynamics interpolated to the coorbital grid, coorbital modes, closed-form Wigner-D rotation `h_inertial[ℓ,m′] = Σ_m D^ℓ_{m′m}(q) h[ℓ,m]`, spline resampling; `mode_sum` with `₋₂Y_ℓm = √((2ℓ+1)/4π) d^ℓ_{m,2}(θ) e^{imφ}`, evaluated at `(inclination, π/2 − phi_ref)`) and the SI generator (`PrecessingGenerator`, `hphc_td`). ForwardDiff and Mooncake pass `integration/GWSurrogates` (`just test-integration-surrogates`); see the package README for the AD rules found while porting (no branches inside recurrences, `element_type` for mixed tuples). It is time-domain only and has no dependency on the GPL packages or `AbstractGWWaveforms` (whose interface is frequency-domain). Model data (Zenodo 10.5281/zenodo.22257361, CC-BY-4.0) is downloaded to a Scratch cache on first use (or read from `GWSURROGATES_DATA_DIR`) and kept in flat integer-indexed vectors, not `Dict`s or closures. Fixtures come from `python/generate_surrogates.py` (own uv project `python/surrogates`, `gwsurrogate==1.2.0`) into `reference/surrogates/`; tests need both the data and the fixtures (CI generates them).
- `GWDistributions`: Distributions.jl-based population models for masses, spins, orientations, delays, and redshifts. Cosmology comes from the external, unregistered `BackgroundCosmology.jl` dependency via a `[sources]` URL.

The waveform packages depend on `GWParameterConversion` through relative `[sources]` paths.
Shared bilby fixture tooling lives in `python/`; generated fixtures go to `reference/` (gitignored).
Waveform scope remains FD only, uniform `AbstractRange` grids starting at 0, detector-frame masses, no tides.
The waveform backends, conversion, detector, series, noise, and matched-filtering packages use GPL-2.0-or-later; `AbstractGWWaveforms`, `GWDistributions` and `GWSurrogates` are MIT.
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
just test-integration-waveforms               # ForwardDiff + Mooncake
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

`packages/GWWaveforms/src/` is generic over `Real`, with AD support restricted to ForwardDiff and Mooncake. Default tests check ForwardDiff against finite differences; `integration/GWWaveforms` checks both backends. Keep kernels type stable and free of `try`/`catch`. Validation runs separately through `validate_parameters` on the host. Grids and generators are fixed configuration; parameter-dependent branches and ordinary indexing are allowed. `maximum_frequency` must be positive. Preserve exact frequency bounds and LAL operation order.

`taylorf2.jl` ports `XLALSimInspiralPNPhasing_F2` and `XLALSimInspiralTaylorF2`. `phenomd.jl` ports `LALSimIMRPhenomD*.c`, with calibrated fits in `phenomd_fits.jl` and QNM tables in `phenomd_qnm.jl`. Its spline uses ordinary segment indexing and its piecewise amplitude/phase use lazy branches. Frequency evaluation uses loops over configured slices, with PhenomD's exclusive upper edge `i < floor(min(0.2 / M, f_max) / Δf)`. Both models reproduce `LALWaveforms` to about 1e-12 relative (largest sample), bar libm differences on Linux. LAL wins over Ripple or papers wherever they differ.

`phenomxas.jl` implements the default IMRPhenomXAS calibration, using the pinned Ripple fits in `phenomxas_fits.jl` and remnant/PN expressions in `phenomxas_utils.jl`, checked against LAL 7.26. Tiny collocation systems use StaticArrays and pivoted solves; differentiate both matrix and RHS. Its cutoff is inclusive at `0.3 / M_sec`; lazy branches evaluate only the active amplitude and phase regions. Direct component-mass differences and a regular near-zero-spin ISCO evaluation avoid endpoint derivative singularities. Preserve these numerical regularizations. `integration/PhenomXASSolve` preserves the earlier solve investigation and is not a runtime dependency.

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
changes also test `GWWaveforms`, and `AbstractGWWaveforms` changes test both waveform packages and `GWSurrogates`. Python tooling changes test all three conversion/waveform
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
