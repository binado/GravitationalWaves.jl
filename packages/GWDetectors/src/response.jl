function _antenna(ra, dec, psi, gmst, D)
    sg, cg = sincos(gmst - ra)
    sd, cd = sincos(dec)
    sp, cp = sincos(psi)
    # Q maps the zero-polarization basis to Earth-fixed coordinates; P rotates that basis by psi.
    Q = @SMatrix [-sg -cg * sd -cg * cd; -cg sg * sd sg * cd; 0 cd -sd]
    P = @SMatrix [cp -sp 0; sp cp 0; 0 0 1]
    R = Q * P
    # Express the detector tensor in the rotated polarization basis.
    response = R' * D * R
    return (; plus = response[1, 1] - response[2, 2], cross = response[1, 2] + response[2, 1])
end

"""
    antenna_pattern(params::NamedTuple, detector)
    antenna_pattern(ra::Number, dec::Number, psi::Number, gmst::Number, detector)

Return `(; plus, cross)` in LAL's polarization convention using a fixed,
long-wavelength detector response. `ra` is right ascension (radians, conventionally
in `[0, 2π)`), `dec` is declination (radians, in `[-π/2, π/2]`), and `psi` is the
polarization angle (radians, periodic modulo `π`). `gmst` is the Greenwich mean
sidereal angle in radians and may be unwrapped. The `params` named tuple supplies `ra`,
`dec`, `psi`, and `geocent_time` in GPS seconds, including fractions; its epoch
follows the conventions of [`greenwich_mean_sidereal_time`](@ref).
Inputs are assumed finite and real; argument values are not validated.
"""
function antenna_pattern(ra::Number, dec::Number, psi::Number, gmst::Number, d::Detector)
    return _antenna(ra, dec, psi, gmst, detector_tensor(d))
end

function antenna_pattern(params::NamedTuple, d::Detector)
    gmst = greenwich_mean_sidereal_time(params.geocent_time)
    return antenna_pattern(params.ra, params.dec, params.psi, gmst, d)
end

function _delay(ra, dec, gmst, vertex)
    sg, cg = sincos(gmst - ra)
    sd, cd = sincos(dec)
    # Reactant's dot overlay does not accept SVector{<:TracedRNumber}.
    return -sum(vertex .* SVector(cd * cg, -cd * sg, sd)) / C_SI
end

"""
    time_delay_from_geocenter(params::NamedTuple, detector)
    time_delay_from_geocenter(ra::Number, dec::Number, gmst::Number, detector)

Arrival time at the detector minus arrival time at Earth's center, in seconds:
`-dot(vertex, direction_to_source)/c`. `ra` is right ascension in radians
(conventionally `[0, 2π)`), `dec` is declination in radians (`[-π/2, π/2]`), and
`gmst` is the Greenwich mean sidereal angle in radians, optionally unwrapped.
`params` supplies `ra`, `dec`, and `geocent_time` in GPS seconds, including fractions;
its epoch follows [`greenwich_mean_sidereal_time`](@ref). Inputs are assumed finite
and real; argument values are not validated.
"""
function time_delay_from_geocenter(ra::Number, dec::Number, gmst::Number, d::Detector)
    return _delay(ra, dec, gmst, d.vertex)
end

function time_delay_from_geocenter(params::NamedTuple, d::Detector)
    gmst = greenwich_mean_sidereal_time(params.geocent_time)
    return time_delay_from_geocenter(params.ra, params.dec, gmst, d)
end

"""
    time_delay(params::NamedTuple, detector1, detector2)

Arrival time at detector 1 minus arrival time at detector 2, in seconds.
`params` supplies right ascension `ra` (radians, conventionally `[0, 2π)`), declination
`dec` (radians, `[-π/2, π/2]`), and `geocent_time` (GPS seconds, including fractions,
following [`greenwich_mean_sidereal_time`](@ref)). Inputs are assumed finite and
real; argument values are not validated.
"""
function time_delay(params::NamedTuple, d1::Detector, d2::Detector)
    gmst = greenwich_mean_sidereal_time(params.geocent_time)
    baseline = d1.vertex .- d2.vertex
    return _delay(params.ra, params.dec, gmst, baseline)
end

"""
    light_travel_time(detector1, detector2)

Euclidean baseline length divided by the speed of light, in seconds.
"""
function light_travel_time(d1::Detector, d2::Detector)
    baseline = d1.vertex .- d2.vertex
    return sqrt(sum(baseline .* baseline)) / C_SI
end

"""
    project_fd(polarizations, sky, detector; epoch, channel=:strain, reference_time=sky.geocent_time)

Project a FrequencySeries with (:plus, :cross) channels into one detector channel.
The destination epoch is explicit; source waveform times are peak-relative.
Raw polarization vectors use a NamedTuple and an additional `frequencies` range
keyword, returning a vector without metadata checks for compiled numerical work.
"""
function project_fd(
        polarizations::NamedTuple, sky, d::Detector;
        frequencies::AbstractRange, epoch::Real, reference_time = sky.geocent_time,
    )
    return _project_fd(first(frequencies), step(frequencies), length(frequencies), epoch, polarizations, sky, d, reference_time)
end

"""Project peak-relative polarization series to a detector series at an explicit destination epoch."""
function project_fd(
        polarizations::FrequencySeries, sky, d::Detector;
        epoch::Real, channel = :strain, reference_time = sky.geocent_time,
    )
    channels(polarizations) == (:plus, :cross) || throw(ArgumentError("polarizations must have channels (:plus, :cross)"))
    # The source FFT is relative to its first sample, which may precede the peak.
    # Add the source offset after subtracting GPS epochs to retain timing precision.
    source = (; plus = view(parent(polarizations), :, 1), cross = view(parent(polarizations), :, 2))
    f = frequencies(polarizations)
    a = _project_fd(first(f), step(f), length(f), epoch, source, sky, d, reference_time, GWSeries.epoch(polarizations))
    return FrequencySeries(a; epoch, frequency_spacing = step(f), channels = (channel,))
end
function _project_fd(f0, df, n, start_time, polarizations, sky, d, reference_time, source_epoch = 0.0)
    hp, hc = polarizations.plus, polarizations.cross
    gmst = greenwich_mean_sidereal_time(reference_time)
    factors = antenna_pattern(sky.ra, sky.dec, sky.psi, gmst, d)
    delay = time_delay_from_geocenter(sky, d)
    dt = (sky.geocent_time - start_time) + delay
    dt += source_epoch
    phase_scale = ((-2 * Float64(π)) * dt) * im
    indices = Base.OneTo(n)
    # One scalar broadcast function prevents Base's range arithmetic shortcuts
    # from constructing a floating-point range inside differentiated code.
    return project_bin.(hp, hc, indices, factors.plus, factors.cross, phase_scale, f0, df)
end

@inline project_bin(hp, hc, i, fp, fc, phase_scale, f0, df) =
    (fp * hp + fc * hc) * exp(phase_scale * (f0 + df * (i - 1)))
