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

The default suites include ForwardDiff checks. Heavy acceptance suites (AD, SNR, live references)
run separately, locally, rather than in CI:

```sh
just test-integration              # all suites, including live LAL and Python references
just test-integration-waveforms    # ForwardDiff and Mooncake waveform gradients
just test-integration-detectors    # detector gradients (ForwardDiff and Mooncake)
```

These recipes require Julia 1.12 or newer and `just`. Set `JULIA='julia +1.12'`
to select a Julia version. The Julia acceptance recipes need no LAL or Python;
the aggregate also runs the two live reference recipes, which require their configured
libraries and Python interpreter. Enzyme and Reactant support is deferred.
Run the affected acceptance suite before submitting numerical or AD changes.

Every package uses independently discoverable TestItems, with TestItemRunner under
`Pkg.test()`. The optional single positional test argument selects items by tags:

| `test_args` | Selection |
|---|---|
| `[]` or `["default"]` | Items without `:integration` (the CI suite; local only for SurrogatePyCallWaveforms) |
| `["integration"]` | Items with `:integration`, excluding `:reference` |
| `["reference"]` | Live items with both `:integration` and `:reference`; fixture-based items tagged only `:reference` run by default |
| `["all"]` | Every item in the package |

Unknown selections and multiple arguments fail before any item executes. For example:

```sh
julia --startup-file=no --check-bounds=yes --project=packages/GWSeries -e 'using Pkg; Pkg.test(test_args=["integration"])'
```

Acceptance tests live under each owning package's `test/` directory. Test dependencies,
including Mooncake where used, are resolved and may be precompiled even for the default
selection, increasing initial installation cost. Numerical packages gain no runtime
dependencies from this migration. PyCall remains confined to SurrogatePyCallWaveforms,
which stays outside workspace membership and CI. Fixture-based reference and lightweight
AD tests remain in the default selection; live reference selections fail if required
libraries, Python imports, or data are unavailable.

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

### Acceptance ownership

Shared cases come from `packages/AbstractWaveformsTests/data/parameters.csv`:
solar masses, Mpc, radians, and dimensionless Cartesian spins in the LAL source frame
at the reference epoch, with inclination measured from orbital angular momentum.
FD generators require aligned magnitude/tilt spins. Surrogates use total mass,
q = m1/m2 (the inverse of bilby's mass ratio), Cartesian spins, and observer azimuth
π/2 - phi_ref. Case selection follows generator capabilities; precessing cases
are never projected onto aligned models.

- [GWSeries](packages/GWSeries/README.md): physical Fourier AD and metadata regression.
- [GWDetectors](packages/GWDetectors/README.md): geometry, timing and projection AD.
- [GWMatchedFiltering](packages/GWMatchedFiltering/README.md): composed likelihoods, vectors and windows.
- [FrequencyDomainWaveforms](packages/FrequencyDomainWaveforms/README.md): waveform AD, SNR, Jacobians and live LAL comparison.
- [SurrogateWaveforms](packages/SurrogateWaveforms/README.md): native kernels, dynamics, waveform AD and SNR.
- [AbstractWaveforms](packages/AbstractWaveforms/README.md): shared harmonic assembly AD.
- [SurrogatePyCallWaveforms](packages/SurrogatePyCallWaveforms/README.md): native-vs-gwsurrogate reference comparison.

### Acceptance SNR conventions

The FD grid is 0:0.5:1024 Hz; generators use 20–1024 Hz. The TD grid has 432 samples
starting at -0.2 seconds with dt=1/2048 s. The reference is the first dynamics node
(`f_ref=0`), and spins are specified in that reference frame. The window has start
and end Planck widths of 0.02 and 0.005 seconds; 80 zero samples are appended after
windowing, for N=512, df=4 Hz. There is no taper-power renormalization. The Fourier
convention is dt*rfft along rows, with no epoch phase factor. An overall epoch
phase would leave this objective unchanged. Native FD outputs are not tapered.

Fixed response coefficients (Fplus,Fcross)=(0.7,0.3) combine both polarizations.
Each model/case is checked with both bundled ET-D (`ET_000A_18_ETD_hlf`) and CE
reference (`cosmic_explorer_strain`) curves. These are one-sided PSDs in strain²/Hz,
read through GWNoise with log-log interpolation. The analysis band is 20–512 Hz
intersected with table/grid coverage. Only active bins are interpolated; no PSD
extrapolation is permitted. GWMatchedFiltering supplies 4df normalization and
excludes DC; the TD metric explicitly receives N and excludes even-length Nyquist.
Native FD grids have no inferred time length.

Configuration, windows, padding buffers, model data, grids, and PSD metrics are
prepared outside differentiation. Only waveform parameters are differentiated.
Mass/distance coordinates are scaled by their fixture values, spins/angles by one.
FD tilt signs are fixed at preparation, so the local signed-spin coordinate is
smooth at zero. Shared PSD preparation lives at
`packages/AbstractWaveformsTests/test/support/snr.jl`. Gradients are compared
componentwise after division by base SNR,
with rtol 1e-8 (FD) / 1e-6 (TD) and atol 1e-8. Values use rtol 1e-10. Finite-difference
spot checks use smooth aligned FD and generic-precession TD cases; boundary and
zero-spin cases compare backend conventions without asserting differentiability
where the physical parameterization has a kink. Every case tests cache reuse at
a nearby distance and the analytic inverse-distance derivative. A phase-sensitive
overlap checks the phase derivative that an SNR norm can miss.

### PSD provenance follow-up

Reuse the maintainer-supplied GWNoise tables verbatim; do not duplicate them here.
Their exact external source/version and licensing attribution were not supplied.
Resolving that attribution remains an explicit follow-up, not a completed claim
of this test migration. GWNoise documents units and validates frequency coverage.

## Formatting

Code is formatted with [Runic](https://github.com/fredrikekre/Runic.jl):

```sh
julia -e 'using Pkg; Pkg.Apps.add("Runic")'
runic --inplace .
```

## Licenses

The root [LICENSE](LICENSE) and every package except `GWSeries`, `AbstractWaveforms`, `GWDistributions`, `GWMatchedFiltering`, `SurrogateWaveforms`, `AbstractWaveformsTests` and `SurrogatePyCallWaveforms` use GPL-2.0-or-later.
`GWSeries` ([LICENSE](packages/GWSeries/LICENSE)), `AbstractWaveforms` ([LICENSE](packages/AbstractWaveforms/LICENSE)), `GWDistributions` ([LICENSE](packages/GWDistributions/LICENSE)), `GWMatchedFiltering` ([LICENSE](packages/GWMatchedFiltering/LICENSE)), `SurrogateWaveforms` ([LICENSE](packages/SurrogateWaveforms/LICENSE)), `AbstractWaveformsTests` ([LICENSE](packages/AbstractWaveformsTests/LICENSE)) and `SurrogatePyCallWaveforms` ([LICENSE](packages/SurrogatePyCallWaveforms/LICENSE)) are MIT.
