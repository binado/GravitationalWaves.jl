# GWDetectors

Native Julia geometry, antenna factors, arrival delays, and frequency-domain projection
for ground-based interferometers. Requires Julia 1.12 and StaticArrays.jl.
The response uses the long-wavelength approximation at one reference time.

```julia
using GWDetectors

detector = Detector("H1")                 # also L1, V1, K1
sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
factors = antenna_pattern(sky, detector)    # (; plus, cross)
delay = time_delay_from_geocenter(sky, detector)  # detector minus geocenter, seconds
baseline_delay = time_delay(sky, detector, Detector("L1")) # H1 minus L1
maximum_delay = light_travel_time(detector, Detector("L1")) # seconds
```

Angles are radians, distances metres, times GPS seconds. Right ascension and
polarization angles are periodic; declination and latitude lie in `[-π/2, π/2]`.
Geocentric delay is `-dot(vertex, direction_to_source)/c`. GMST is unwrapped and
uses LAL's J2000 convention, with separate whole and fractional seconds and a
pinned leap table through 2017-01-01. Inputs in `[0, 2147483647]` are supported;
after the last table entry GPS−UTC is assumed to remain 18 seconds. This follows
LAL's UTC-based approximation without a UT1 correction. At a positive leap second,
23:59:60 and the following 00:00:00 map to the same Julian day.

## Composition with waveform generators

Both waveform backends return the `(; plus, cross)` accepted by `project_fd`:

```julia
using GWWaveforms, LALWaveforms, GWDetectors

f = range(0, 512; step = 0.25)
p = (; mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
      theta_jn = 0.4, phase = 1.3)
for generator in (IMRPhenomD(), LALWaveformGenerator("IMRPhenomD"))
    polarizations = hphc_fd(f, p, generator)
    strain = project_fd(f, polarizations, sky, detector;
                        start_time = 1126259460.0)
end
```

`project_fd` allocates a complex vector in input order for any finite nonnegative
frequency vector (ranges, irregular, unsorted, repeated, and empty grids included).
Both polarization vectors must have the same axes as the frequencies. It applies
no additional frequency mask. The expression is

```math
(F_+ h_+ + F_\times h_\times)
\exp\{-2\pi i f [(t_c-t_\mathrm{start})+\tau(t_c)]\}.
```

The epoch difference is evaluated before adding the geometric delay. Antenna
factors use `reference_time` (defaults to `sky.geocent_time`), while the delay
always uses `sky.geocent_time`. Supply `reference_time` explicitly to hold the
response fixed across draws. Delay methods need only `ra`, `dec`, `geocent_time`.
For fixed GMST, use `antenna_pattern(ra, dec, psi, gmst, detector)` and
`time_delay_from_geocenter(ra, dec, gmst, detector)`.

## Custom geometry

```julia
custom = Detector("custom"; latitude = 0.6, longitude = -1.1, elevation = 50.0,
    xarm_azimuth = 0.2, yarm_azimuth = 0.2 + π/3,
    xarm_tilt = 0.001, yarm_tilt = -0.002, arm_lengths = (3000.0, 3000.0))
cartesian = Detector("copy", custom.vertex, custom.xarm, custom.yarm, custom.arm_lengths)
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
in their docstrings and do not validate values. Projection checks polarization
vector types and matching axes. Geometry constructors accept `check=false` to
bypass numerical geometry validation inside compiled calls with prevalidated
inputs; vector sizes and numeric element types are still checked. GPS derivatives
are tested away from integer-second rounding
boundaries and leap transitions. Scalar geometry uses StaticArrays; projection uses
broadcast so device arrays stay in compiled execution.

```julia
using Reactant, Enzyme
Reactant.set_default_backend("cpu")

function projected(angles, f, hp, hc)
    Reactant.@allowscalar begin
        source = (; ra = angles[1], dec = angles[2], psi = angles[3],
                    geocent_time = angles[4])
        project_fd(f, (; plus = hp, cross = hc), source, Detector("H1");
                   start_time = 1126259460.0)
    end
end
angles = Reactant.to_rarray([sky.ra, sky.dec, sky.psi, sky.geocent_time])
frequencies = Reactant.to_rarray([0.0, 20.0, 35.5])
hp = Reactant.to_rarray(ComplexF64[0, 1, 0.5im])
hc = Reactant.to_rarray(ComplexF64[0, 0.3im, 0.2])
executable = Reactant.@compile projected(angles, frequencies, hp, hc)
strain = executable(angles, frequencies, hp, hc)
```

`just test-integration-detectors` exercises reused CPU executables with changing
sky, GPS, Cartesian/geodetic geometry, and polarization/frequency arrays, including
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
