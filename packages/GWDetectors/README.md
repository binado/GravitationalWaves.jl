# GWDetectors

Native Julia geometry, antenna factors, arrival delays, and frequency-domain projection
for ground-based interferometers. Requires Julia 1.12, StaticArrays.jl, and GWSeries.
The response uses the long-wavelength approximation at one reference time.

```julia
using GWDetectors

detector = detectors.H1                   # also detectors.L1, detectors.V1, detectors.K1
sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
factors = antenna_pattern(sky, detector)    # (; plus, cross)
delay = time_delay_from_geocenter(sky, detector)  # detector minus geocenter, seconds
baseline_delay = time_delay(sky, detector, detectors.L1) # H1 minus L1
maximum_delay = light_travel_time(detector, detectors.L1) # seconds
```

`detectors` is an immutable named tuple of predefined `Detector` instances;
use `detectors.H1` or `detectors[:H1]`. Construct custom geometry with `Detector`
as shown below.

Angles are radians, distances metres, times GPS seconds. Right ascension and
polarization angles are periodic; declination and latitude lie in `[-π/2, π/2]`.
Geocentric delay is `-dot(vertex, direction_to_source)/c`. GMST is unwrapped and
uses LAL's J2000 convention, with separate whole and fractional seconds and a
pinned leap table through 2017-01-01. Inputs in `[0, 2147483647]` are supported;
after the last table entry GPS−UTC is assumed to remain 18 seconds. This follows
LAL's UTC-based approximation without a UT1 correction. At a positive leap second,
23:59:60 and the following 00:00:00 map to the same Julian day.

## Sensitivities and detector designs

`detectors` also includes triangular ET sites `E1`–`E3`, ET L sites
`S1`, `S2`, `R1`, `R2`, and CE sites `C1`, `C2`. Their inline definitions are transcribed from the supplied
[geometry catalog](data/geometry.toml), with nominal full arm
lengths of 10 km, 15 km, and 40 km respectively. `S1`/`R1` are aligned;
`S2`/`R2` are misaligned. `C1`/`C2` use the Hanford/Livingston sites.
The original four detectors retain their pinned LAL Cartesian constants.

GWNoise bundles twelve text tables of one-sided PSDs
in strain²/Hz, with frequency in Hz. Load them explicitly with GWNoise:

```julia
using GWDetectors, GWNoise

curve = read_psd(psd_path(sensitivities.H1))
value = psd(100.0, curve)
curve_names = available_psds()  # all designs, including alternatives
alternative = read_psd(psd_path(:ET_COBA_20km); interpolation = :loglog)
```

`psd_path` accepts symbols or strings. Detector defaults in `sensitivities` come from the supplied
[sensitivity catalog](data/sensitivity.toml): H1/L1 use A+,
V1 uses Virgo O5 low, K1 uses KAGRA 128 Mpc, E1–E3 use ET COBA 10 km,
S1/S2/R1/R2 use ET COBA 15 km, and C1/C2 use the CE reference curve.
Curve names identify sensitivity designs independently of geometry; choosing
an alternative PSD does not change arm lengths. GWNoise controls interpolation
and frequency coverage; its default is to throw outside the tabulated band.
GWNoise provides local paths; GWDetectors maps detector names to curve names
without depending on GWNoise at runtime.
The supplied catalog comments retain their original project context.

## Composition with waveform generators

Use a series prototype to generate FrequencySeries polarizations with channels `(:plus, :cross)`:

```julia
using FrequencyDomainWaveforms, GWDetectors, GWSeries
f = 0.0:0.25:512.0
p = (; mass_1=36.0, mass_2=29.0, luminosity_distance=500.0, theta_jn=0.4, phase=1.3)
sky = (; ra=1.2, dec=-0.4, psi=0.7, geocent_time=1126259462.25)
prototype = FrequencySeries(zeros(ComplexF64, length(f)); epoch=1126259460.0, frequency_spacing=step(f), channels=(:H1,))
polarizations = hphc_fd(prototype, p, IMRPhenomD())
strain = project_fd(polarizations, sky, detectors.H1; epoch=1126259460.0, channel=:H1)
```

The destination epoch is explicit. Output channels default to `(:strain,)`; supply
`channel` for a detector label. The source epoch is peak-relative and is included
in the projection phase, so FFT-converted TD waveforms retain their first-sample
reference. Epoch subtraction precedes source-offset addition to retain GPS precision.

For unchecked compiled numerical work, pass raw polarization vectors as
`project_fd((;plus=hp,cross=hc), sky, detector; frequencies=f, epoch=start_time,
reference_time=sky.geocent_time)`. This returns an ordinary vector. The phase kernel
broadcasts over integer frequency-bin indices and preserves the existing LAL/bilby
operation order. `reference_time` sets antenna factors; arrival delays use
`sky.geocent_time`. This is the fixed long-wavelength response.

## Custom geometry

```julia
custom = Detector(; latitude = 0.6, longitude = -1.1, elevation = 50.0,
    xarm_azimuth = 0.2, yarm_azimuth = 0.2 + π/3,
    xarm_tilt = 0.001, yarm_tilt = -0.002, arm_lengths = (3000.0, 3000.0))
cartesian = Detector(custom.vertex, custom.xarm, custom.yarm, custom.arm_lengths)
D = detector_tensor(custom)  # 3×3 SMatrix
```

Geodetic construction uses WGS-84 (radius 6378137 m, inverse flattening
298.257223563). Azimuth increases clockwise from north; tilt is positive above
the local horizontal. Arms may be nonorthogonal but must be nonparallel unit
vectors. Unit-norm validation allows absolute error `1e-6` to preserve KAGRA's
published precision. Cartesian directions are never renormalized. Arm lengths
are full lengths, stored as metadata; they do not enter the long-wavelength
response tensor `(xarm ⊗ xarm - yarm ⊗ yarm)/2`. Input types are promoted without
forcing Float64. Geometry is real-valued; `Detector{T<:Number}` also accommodates
Reactant scalar tracing. Cartesian vectors and arm lengths are stored as
`SVector`s; `detector_tensor` returns an `SMatrix`, indexed as `D[i, j]`.

The antenna response uses `A = R' * D * R`, where `R` rotates the polarization
basis into Earth-fixed coordinates, including the polarization angle `psi`.
The factors are `Fplus = A[1, 1] - A[2, 2]` and
`Fcross = A[1, 2] + A[2, 1]`, in LAL's convention.

## Differentiation and Reactant CPU

ForwardDiff, Enzyme, Mooncake, and compiled Enzyme through Reactant are tested.
Response, timing, and projection functions assume the argument domains documented
in their docstrings and do not validate values. Projection assumes prevalidated polarization vectors and a fixed domain. Geometry constructors accept `check=false` to
bypass numerical geometry validation inside compiled calls with prevalidated
inputs; vector sizes and numeric element types are still checked. GPS derivatives
are tested away from integer-second rounding
boundaries and leap transitions. Scalar geometry uses StaticArrays; projection uses
broadcast so device arrays stay in compiled execution.

```julia
using Reactant, Enzyme, GWSeries
Reactant.set_default_backend("cpu")

const PROJECTION_FREQUENCIES = 0.0:20.0:40.0
const PROJECTION_EPOCH = 1126259460.0
function projected(angles, hp, hc)
    Reactant.@allowscalar begin
        source = (;ra=angles[1],dec=angles[2],psi=angles[3],geocent_time=angles[4])
        project_fd((;plus=hp,cross=hc),source,detectors.H1; frequencies=PROJECTION_FREQUENCIES, epoch=PROJECTION_EPOCH)
    end
end
angles = Reactant.to_rarray([sky.ra, sky.dec, sky.psi, sky.geocent_time])
hp_host = ComplexF64[0,1,0.5im]
hc_host = ComplexF64[0,0.3im,0.2]
hp, hc = Reactant.to_rarray.((hp_host,hc_host))
executable = Reactant.@compile projected(angles,hp,hc)
strain = executable(angles,hp,hc)
```

`just test-integration-detectors` exercises reused CPU executables with changing
sky, GPS, Cartesian/geodetic geometry, and polarization arrays at fixed frequency domains, including
phase-sensitive compiled gradients. GPU compilation and compilation of waveform
generators are outside this package's acceptance scope.

## References and license

GPL-2.0-or-later; see [LICENSE](LICENSE). Antenna, timing, and sidereal formulas
are ported from LAL's `DetResponse.c`, `TimeDelay.c`, and `XLALSiderealTime.c`
(Copyright Jolien Creighton, Kipp Cannon, David Chin, Steven Fairhurst, Alexander
Dietz, Drew Keppel, and other LAL contributors). The pinned built-in constants
are the `LHO_4K`, `LLO_4K`, `VIRGO`, and `KAGRA` macros in `LALDetectors.h` from
LAL 7.7.1 / lalsuite 7.26.16. Full lengths are twice the published arm midpoints.
Sources: [LAL detector constants](https://lscsoft.docs.ligo.org/lalsuite/dev/lal/group___detector_constants.html),
[LAL source](https://git.ligo.org/lscsoft/lalsuite), and
[Anderson et al., PRD 63 042003](https://doi.org/10.1103/PhysRevD.63.042003).

Generate independent fixtures with
`uv run --project python python python/generate_detectors.py`. The generator
records package versions, Cartesian inputs, native double tensors, rounded
REAL4 tensors, fixed-GMST responses, leap-transition GMST values, and complete
Bilby projections with matched Cartesian geometry and reference epochs. Missing
fixtures warn locally and fail when `CI=true`.

Finite-arm transfer functions, rotating time-domain response, extra polarizations,
noise, calibration, and a network container remain future work.

## Series migration

Projection accepts polarization series and returns detector series. For repeated AD
work, retain fixed frequency/epoch keywords and plain sample arrays. Metadata is
configuration; runtime-changing grids and FFT differentiation are outside acceptance.
