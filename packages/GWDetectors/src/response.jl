function _check_sky(ra, dec, angles...)
    all(_finite_real, (ra, dec, angles...)) || throw(ArgumentError("sky angles must be finite and real"))
    abs(dec) <= π / 2 || throw(ArgumentError("declination must be in [-π/2, π/2]"))
    return nothing
end

# XLALComputeDetAMResponse: keep both contractions in one traversal.
function _antenna(ra, dec, psi, gmst, D)
    sg, cg = sincos(gmst - ra)
    sd, cd = sincos(dec)
    sp, cp = sincos(psi)
    X = (-cp * sg - sp * cg * sd, -cp * cg + sp * sg * sd, sp * cd)
    Y = (sp * sg - cp * cg * sd, sp * cg + cp * sg * sd, cp * cd)
    plus = zero(X[1] * D[1][1])
    cross = plus
    for i in 1:3
        DX = _dot(D[i], X)
        DY = _dot(D[i], Y)
        plus += X[i] * DX - Y[i] * DY
        cross += X[i] * DY + Y[i] * DX
    end
    return (; plus, cross)
end

"""
    antenna_pattern(sky, detector; check=true)
    antenna_pattern(ra, dec, psi, gmst, detector; check=true)

Return `(; plus, cross)` in LAL's polarization convention. `sky` supplies `ra`,
`dec`, `psi` (radians), and `geocent_time` (GPS seconds). The explicit form takes
GMST in radians. Both use a fixed, long-wavelength detector response.
"""
function antenna_pattern(ra, dec, psi, gmst, d::Detector; check = true)
    check && _check_sky(ra, dec, psi, gmst)
    return _antenna(ra, dec, psi, gmst, detector_tensor(d))
end

function antenna_pattern(sky, d::Detector; check = true)
    gmst = greenwich_mean_sidereal_time(sky.geocent_time; check)
    return antenna_pattern(sky.ra, sky.dec, sky.psi, gmst, d; check)
end

function _delay(ra, dec, gmst, vertex)
    sg, cg = sincos(gmst - ra)
    sd, cd = sincos(dec)
    return -_dot(vertex, (cd * cg, -cd * sg, sd)) / C_SI
end

"""
    time_delay_from_geocenter(sky, detector; check=true)
    time_delay_from_geocenter(ra, dec, gmst, detector; check=true)

Arrival time at the detector minus arrival time at Earth's center, in seconds:
`-dot(vertex, direction_to_source)/c`. `sky` needs only `ra`, `dec`, `geocent_time`.
"""
function time_delay_from_geocenter(ra, dec, gmst, d::Detector; check = true)
    check && _check_sky(ra, dec, gmst)
    return _delay(ra, dec, gmst, d.vertex)
end

function time_delay_from_geocenter(sky, d::Detector; check = true)
    gmst = greenwich_mean_sidereal_time(sky.geocent_time; check)
    return time_delay_from_geocenter(sky.ra, sky.dec, gmst, d; check)
end

"""
    time_delay(sky, detector1, detector2; check=true)

Arrival time at detector 1 minus arrival time at detector 2, in seconds.
"""
function time_delay(sky, d1::Detector, d2::Detector; check = true)
    gmst = greenwich_mean_sidereal_time(sky.geocent_time; check)
    check && _check_sky(sky.ra, sky.dec)
    baseline = ntuple(i -> d1.vertex[i] - d2.vertex[i], 3)
    return _delay(sky.ra, sky.dec, gmst, baseline)
end

"""
    light_travel_time(detector1, detector2)

Euclidean baseline length divided by the speed of light, in seconds.
"""
function light_travel_time(d1::Detector, d2::Detector)
    baseline = ntuple(i -> d1.vertex[i] - d2.vertex[i], 3)
    return sqrt(_dot(baseline, baseline)) / C_SI
end

"""
    project_fd(f, polarizations, sky, detector;
               start_time, reference_time=sky.geocent_time, check=true)

Allocate `(Fplus*hplus + Fcross*hcross) * exp(-2π*im*f*dt)`, with
`dt = (sky.geocent_time - start_time) + time_delay_from_geocenter(sky, detector)`.
Antenna factors are evaluated at `reference_time`; delay at `sky.geocent_time`.
Epochs are subtracted before adding the delay to preserve timing precision.

Accept any finite nonnegative frequency vector, including ranges, irregular or
empty grids, in input order. Both polarization vectors must have the same axes as
`f`; no frequency mask is applied. `sky` provides `ra`, `dec`, `psi`, `geocent_time`.
`check=false` skips value validation for prevalidated compiled inputs, retaining
vector and axes checks. Broadcast preserves compiled array execution.
"""
function project_fd(
        f::AbstractVector, polarizations, sky, d::Detector;
        start_time, reference_time = sky.geocent_time, check = true
    )
    hp, hc = polarizations.plus, polarizations.cross
    hp isa AbstractVector && hc isa AbstractVector || throw(ArgumentError("polarizations must be vectors"))
    axes(f) == axes(hp) == axes(hc) || throw(DimensionMismatch("frequency and polarization axes must match"))
    if check
        all(x -> _finite_real(x) && x >= 0, f) || throw(ArgumentError("frequencies must be finite, real and nonnegative"))
        all(isfinite, hp) && all(isfinite, hc) || throw(ArgumentError("polarizations must be finite"))
        _finite_real(start_time) || throw(ArgumentError("start_time must be finite and real"))
    end
    gmst = greenwich_mean_sidereal_time(reference_time; check)
    factors = antenna_pattern(sky.ra, sky.dec, sky.psi, gmst, d; check)
    delay = time_delay_from_geocenter(sky, d; check)
    dt = (sky.geocent_time - start_time) + delay
    phase_scale = ((-2 * Float64(π)) * dt) * im
    return @. (factors.plus * hp + factors.cross * hc) * exp(phase_scale * f)
end
