# GravitationalWaves.jl

A collection of composable Julia packages for gravitational-wave modeling. The repository is a
Julia workspace; each package has its own dependencies, tests, version, and public API.

| Package | Purpose | License |
|---------|---------|---------|
| [AbstractGWWaveforms](packages/AbstractGWWaveforms) | Dependency-free waveform-generator interface shared by the waveform packages | GPL-2.0-or-later |
| [GWWaveforms](packages/GWWaveforms) | Native, differentiable TaylorF2 and IMRPhenomD waveforms | GPL-2.0-or-later |
| [LALWaveforms](packages/LALWaveforms) | Direct calls to LALSimulation frequency-domain waveforms | GPL-2.0-or-later |
| [GWParameterConversion](packages/GWParameterConversion) | Differentiable mass, spin, and frame conversions | GPL-2.0-or-later |
| [GWDetectors](packages/GWDetectors) | Native, differentiable detector geometry, antenna patterns, timing, and FD projection | GPL-2.0-or-later |
| [GWDistributions](packages/GWDistributions) | Population distributions for masses, spins, orientations, and redshifts | MIT |

## Installation

Install individual packages by subdirectory. The packages are currently unregistered, so include
unregistered dependencies explicitly in the active environment:

```julia
using Pkg

repo = "https://github.com/binado/GravitationalWaves.jl"
Pkg.add([
    PackageSpec(url = repo, subdir = "packages/GWParameterConversion"),
    PackageSpec(url = repo, subdir = "packages/GWWaveforms"),
    PackageSpec(url = repo, subdir = "packages/GWDetectors"),
])
```

For `LALWaveforms`, replace `packages/GWWaveforms` above with `packages/LALWaveforms`.
For population distributions:

```julia
Pkg.add([
    PackageSpec(name = "Distributions"),
    PackageSpec(url = "https://github.com/binado/BackgroundCosmology.jl"),
    PackageSpec(url = repo, subdir = "packages/GWDistributions"),
])
```

`BackgroundCosmology` remains an external dependency. There is no umbrella Julia module.

## LALWaveforms

Frequency-domain compact-binary waveforms from [LALSimulation](https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/),
called directly through `ccall`, with [bilby](https://bilby-dev.github.io/bilby/)'s parameter names, units and
parameter-to-LAL conversions.

```julia
using LALWaveforms

generator = LALWaveformGenerator("IMRPhenomXPHM"; reference_frequency = 50.0, minimum_frequency = 20.0)
f = range(0, 1024; step = 0.25)
params = (
    mass_1 = 36.0, mass_2 = 29.0,          # detector-frame solar masses
    luminosity_distance = 500.0,           # Mpc
    theta_jn = 0.4, phase = 1.3,
    a_1 = 0.5, a_2 = 0.3, tilt_1 = 1.0, tilt_2 = 0.6, phi_12 = 0.3, phi_jl = 4.0,
)
(; plus, cross) = hphc_fd(f, params, generator)
```

`LALWaveformGenerator` builds the LAL waveform dictionary once; reuse it across calls. Optional parameters
(`a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, eccentricity`) default to 0, and unknown keys throw. Instead of
`mass_1`/`mass_2` you may pass any other mass parameters that
[GWParameterConversion](packages/GWParameterConversion) can convert (for example `chirp_mass` and `mass_ratio`).

Scope for now: frequency domain only, uniform `AbstractRange` grids starting at 0, detector-frame masses,
no tidal parameters.

## GWWaveforms: native, differentiable waveforms

[`packages/GWWaveforms`](packages/GWWaveforms) provides the same interface in pure Julia (no LAL install),
for TaylorF2 and IMRPhenomD. The code is generic over `Real`, so it works with ForwardDiff, Enzyme and
Mooncake, and its tests check it against `LALWaveforms` to about 1e-12.

```julia
using GWWaveforms

generator = IMRPhenomD()
(; plus, cross) = hphc_fd(range(0, 1024; step = 0.25), params, generator)
```

Both packages implement the interface of [`AbstractGWWaveforms`](packages/AbstractGWWaveforms): a concrete
`AbstractWaveformGenerator` (`TaylorF2()`, `IMRPhenomD()`, `LALWaveformGenerator("...")`) and the one generic
`hphc_fd(f, params, generator)`, plus `validate_parameters(f, params, generator)` for the parameter
checks. `GWWaveforms.hphc_fd` does not check parameter values, so that it can be compiled with Reactant (see its
[README](packages/GWWaveforms)); run `validate_parameters` on the host first. `LALWaveforms.hphc_fd` still validates.
Loading both with `using` is unambiguous, and generic code can take any backend.

## Detector projection

[`GWDetectors`](packages/GWDetectors) composes with both waveform backends:

```julia
using GWDetectors

sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
detector = Detector("H1")  # H1, L1, V1, K1 or custom Cartesian/geodetic geometry
for generator in (IMRPhenomD(), LALWaveformGenerator("IMRPhenomD"))
    polarizations = hphc_fd(f, params, generator)
    strain = project_fd(f, polarizations, sky, detector; start_time = 1126259460.0)
end
```

It uses a long-wavelength response evaluated at one reference time and supports
arbitrary frequency vectors. The package has no runtime dependencies and supports
ForwardDiff, Enzyme, Mooncake, and Reactant CPU compilation. See its
[README](packages/GWDetectors/README.md) for timing conventions, custom tilted or
nonorthogonal arms, compiled Reactant examples, and leap-table coverage.

## Pointing the package at LAL

LAL is not bundled. Provide the directory holding `liblal` and `liblalsimulation` in one of two ways
(the environment variable wins):

- `LALWAVEFORMS_LIBDIR=/path/to/libs`
- `LALWaveforms.set_libdir!("/path/to/libs")`, saved as a preference, or
  `LALWaveforms.set_libdir_from_python!()` to take it from the `lalsuite` wheel in the active Python environment.

The `lalsuite` wheels carry the libraries in `lal/.dylibs/` (macOS) or `lalsuite.libs/` (Linux).
Loading is lazy: the package loads without LAL and the first waveform call explains how to configure it.

## Development and tests

From the repository root, instantiate the workspace and test each package:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(; workspace = true)'
julia --project=packages/AbstractGWWaveforms -e 'using Pkg; Pkg.test()'
julia --project=packages/GWDistributions -e 'using Pkg; Pkg.test()'
julia --project=packages/GWParameterConversion -e 'using Pkg; Pkg.test()'
julia --project=packages/GWWaveforms -e 'using Pkg; Pkg.test()'
julia --project=packages/GWDetectors -e 'using Pkg; Pkg.test()'
```

The waveform tests require the reference setup below for full coverage.

The default suites include ForwardDiff checks. Heavy Reactant/Enzyme conversion tests
and Enzyme/Mooncake waveform and detector tests run separately, locally, rather than in CI:

```sh
just test-integration              # all suites; installs their dependencies on first use
just test-integration-conversion   # Reactant CPU compilation and Enzyme gradients
just test-integration-waveforms    # all three waveform AD backends
just test-integration-detectors    # detector gradients and Reactant CPU acceptance
```

These recipes require Julia 1.12 or newer and `just`, and need no LAL or Python setup.
Set `JULIA='julia +1.12'` to select a Julia version. Their standalone environments live
under `integration/`, outside the root workspace; ordinary `Pkg.test()` does not
install or precompile Enzyme, Reactant, or Mooncake. Run the affected integration
suite before submitting changes to the numerical kernels or AD integration tests.

The reference tests compare against bilby output generated locally (not committed):

```sh
uv sync --project python
uv run --project python python python/generate.py
uv run --project python python python/generate_detectors.py
export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)
julia --project=packages/LALWaveforms -e 'using Pkg; Pkg.test()'
```

Julia loads `libopenlibm` globally, so on Linux LAL's libm calls resolve to it instead of glibc. Output is bit-identical to Python on macOS, but differs at the ulp level on Linux (up to ~1e-6 for IMRPhenomT, whose time-domain pipeline amplifies it), so `cases.json` carries a per-case `rtol`.

To get bit-for-bit agreement with Python on Linux, start Julia with glibc's libm forced ahead of openlibm (`LD_PRELOAD` is read by the dynamic loader at process start, so it must be set before launching Julia; the path varies by distro, and it affects every native library in the process):

```sh
LD_PRELOAD=/lib/x86_64-linux-gnu/libm.so.6 julia --project
```

Without fixtures the reference tests warn and skip; with `CI=true` they fail.

## Formatting

Code is formatted with [Runic](https://github.com/fredrikekre/Runic.jl):

```sh
julia -e 'using Pkg; Pkg.Apps.add("Runic")'
runic --inplace .
```

## Licenses

The root [LICENSE](LICENSE) and the waveform, conversion, and detector packages use GPL-2.0-or-later.
`GWDistributions` retains its [MIT license](packages/GWDistributions/LICENSE).
