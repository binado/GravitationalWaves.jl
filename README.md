# GravitationalWaves.jl

A collection of composable Julia packages for gravitational-wave modeling. The repository is a
Julia workspace; each package has its own dependencies, tests, version, and public API.

| Package | Purpose | License |
|---------|---------|---------|
| [AbstractWaveforms](packages/AbstractWaveforms) | Waveform-generator interface using GWSeries shared by the waveform packages | MIT |
| [FrequencyDomainWaveforms](packages/FrequencyDomainWaveforms) | Native, differentiable TaylorF2, IMRPhenomD and IMRPhenomXAS waveforms | GPL-2.0-or-later |
| [LALWaveforms](packages/LALWaveforms) | Direct calls to LALSimulation frequency-domain waveforms | GPL-2.0-or-later |
| [GWParameterConversion](packages/GWParameterConversion) | Differentiable mass, spin, and frame conversions | GPL-2.0-or-later |
| [GWDetectors](packages/GWDetectors) | Native, differentiable detector geometry, antenna patterns, timing, and FD projection | GPL-2.0-or-later |
| [GWDistributions](packages/GWDistributions) | Population distributions for masses, spins, orientations, and redshifts | MIT |
| [GWSeries](packages/GWSeries) | Sample/channel series, reconstructible metadata, and physical AbstractFFTs plans | MIT |
| [GWNoise](packages/GWNoise) | Read PSD/ASD tables and sample stationary Gaussian strain noise | GPL-2.0-or-later |
| [GWMatchedFiltering](packages/GWMatchedFiltering) | Joint frequency-domain metrics, inner products, and SNRs with correlated detector noise | MIT |
| [SurrogateWaveforms](packages/SurrogateWaveforms) | Differentiable NRSur7dq4 numerical-relativity surrogate (in progress) | MIT |
| [AbstractWaveformsTests](packages/AbstractWaveformsTests) | Shared contract, reference-generator and gradient test suites for AbstractWaveforms backends | MIT |
| [SurrogatePyCallWaveforms](packages/SurrogatePyCallWaveforms) | PyCall adapter around upstream `gwsurrogate`, used as a reference generator (local only) | MIT |

## Installation

Install individual packages by subdirectory. The packages are currently unregistered, so include
unregistered dependencies explicitly in the active environment:

```julia
using Pkg

repo = "https://github.com/binado/GravitationalWaves.jl"
Pkg.add([
    PackageSpec(url = repo, subdir = "packages/GWParameterConversion"),
    PackageSpec(url = repo, subdir = "packages/FrequencyDomainWaveforms"),
    PackageSpec(url = repo, subdir = "packages/GWSeries"),
    PackageSpec(url = repo, subdir = "packages/GWDetectors"),
])
```

For `LALWaveforms`, replace `packages/FrequencyDomainWaveforms` above with `packages/LALWaveforms`.
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
h = hphc_fd(f, params, generator)  # ordinary matrix: plus/cross columns
```

`LALWaveformGenerator` builds the LAL waveform dictionary once; reuse it across calls. Optional parameters
(`a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, eccentricity`) default to 0, and unknown keys throw. Instead of
`mass_1`/`mass_2` you may pass any other mass parameters that
[GWParameterConversion](packages/GWParameterConversion) can convert (for example `chirp_mass` and `mass_ratio`).

Scope for now: frequency domain only, uniform `AbstractRange` grids starting at 0, detector-frame masses,
no tidal parameters.

## FrequencyDomainWaveforms: native, differentiable waveforms

[`packages/FrequencyDomainWaveforms`](packages/FrequencyDomainWaveforms) provides the same interface in pure Julia (no LAL install),
for TaylorF2, IMRPhenomD and IMRPhenomXAS. The code is generic over `Real`, with AD support for ForwardDiff and
Mooncake, and its tests check it against `LALWaveforms` to about 1e-12.

```julia
using FrequencyDomainWaveforms

generator = IMRPhenomD()
h = hphc_fd(range(0, 1024; step = 0.25), params, generator)  # plus/cross columns
```

Both packages implement the interface of [`AbstractWaveforms`](packages/AbstractWaveforms): a concrete
`AbstractWaveformGenerator` (`TaylorF2()`, `IMRPhenomD()`, `IMRPhenomXAS()`, `LALWaveformGenerator("...")`) and the one generic
`hphc_fd(f, params, generator)`, plus `validate_parameters(f, params, generator)` for the parameter
checks. `FrequencyDomainWaveforms.hphc_fd` does not check parameter values; run `validate_parameters` on the host first
(see its [README](packages/FrequencyDomainWaveforms)). `LALWaveforms.hphc_fd` still validates.
Loading both with `using` is unambiguous, and generic code can take any backend.

## Detector projection

[`GWDetectors`](packages/GWDetectors) composes with both waveform backends:

```julia
using GWDetectors, GWSeries

f = 0.0:0.25:512.0
sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
detector = detectors.H1  # built-ins: H1, L1, V1, K1; use Detector(...) for custom geometry
for generator in (IMRPhenomD(), LALWaveformGenerator("IMRPhenomD"))
    prototype = FrequencySeries(zeros(ComplexF64, length(f)); epoch=1126259460.0, frequency_spacing=step(f), channels=(:H1,))
    polarizations = hphc_fd(prototype, params, generator)
    strain = project_fd(polarizations, sky, detector; epoch=1126259460.0, channel=:H1)
end
```

It uses a long-wavelength response evaluated at one reference time and supports
uniform frequency domains. The package uses StaticArrays.jl and supports
ForwardDiff and Mooncake. See its
[README](packages/GWDetectors/README.md) for timing conventions, custom tilted or
nonorthogonal arms, and leap-table coverage.

For noise-weighted comparisons of projected strains, see
[GWMatchedFiltering](packages/GWMatchedFiltering/README.md). It works on plain ordered multi-detector
Fourier arrays and either independent PSDs or full Hermitian cross-spectral
matrices. Whitening-factor preparation (`GWNoise.whitening_factor`) is separate from repeated
inner products, SNRs, and fixed-PSD Gaussian likelihood ratios. Layout measurements have a standalone
[benchmark harness](benchmarks/GWMatchedFiltering).

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
julia --project=packages/AbstractWaveforms -e 'using Pkg; Pkg.test()'
julia --project=packages/GWDistributions -e 'using Pkg; Pkg.test()'
julia --project=packages/GWParameterConversion -e 'using Pkg; Pkg.test()'
julia --project=packages/FrequencyDomainWaveforms -e 'using Pkg; Pkg.test()'
julia --project=packages/GWDetectors -e 'using Pkg; Pkg.test()'
```

The waveform tests require the reference setup below for full coverage.

The default suites include ForwardDiff checks. Heavy ForwardDiff/Mooncake waveform, detector and inference
tests run separately, locally, rather than in CI:

```sh
just test-integration              # all suites; installs their dependencies on first use
just test-integration-waveforms    # ForwardDiff and Mooncake waveform gradients
just test-integration-detectors    # detector gradients (ForwardDiff and Mooncake)
```

These recipes require Julia 1.12 or newer and `just`, and need no LAL or Python setup.
Set `JULIA='julia +1.12'` to select a Julia version. Their standalone environments live
under `integration/`, outside the root workspace; ordinary `Pkg.test()` does not
install or precompile Mooncake. Enzyme and Reactant support is deferred. Run the affected integration
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

The root [LICENSE](LICENSE) and every package except `GWSeries`, `AbstractWaveforms`, `GWDistributions`, `GWMatchedFiltering`, `SurrogateWaveforms`, `AbstractWaveformsTests` and `SurrogatePyCallWaveforms` use GPL-2.0-or-later.
`GWSeries` ([LICENSE](packages/GWSeries/LICENSE)), `AbstractWaveforms` ([LICENSE](packages/AbstractWaveforms/LICENSE)), `GWDistributions` ([LICENSE](packages/GWDistributions/LICENSE)), `GWMatchedFiltering` ([LICENSE](packages/GWMatchedFiltering/LICENSE)), `SurrogateWaveforms` ([LICENSE](packages/SurrogateWaveforms/LICENSE)), `AbstractWaveformsTests` ([LICENSE](packages/AbstractWaveformsTests/LICENSE)) and `SurrogatePyCallWaveforms` ([LICENSE](packages/SurrogatePyCallWaveforms/LICENSE)) are MIT.

Waveform numerical acceptance uses a shared physical parameter CSV, reference-generator
comparisons, and optimal-SNR backend comparisons; see [the suite guide](integration/README.md).
