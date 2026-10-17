_finite_real(x) = x isa Real && isfinite(x)

function _check_geometry(vertex, xarm, yarm, lengths)
    all(_finite_real, (vertex..., xarm..., yarm..., lengths...)) ||
        throw(ArgumentError("geometry must be finite and real"))
    # KAGRA's published directions carry only seven decimal places.
    all(u -> isapprox(dot(u, u), 1; atol = 1.0e-6, rtol = 0), (xarm, yarm)) ||
        throw(ArgumentError("arm directions must be unit vectors (tolerance 1e-6)"))
    any(!iszero, cross(xarm, yarm)) || throw(ArgumentError("arm directions must not be parallel"))
    all(>(0), lengths) || throw(ArgumentError("full arm lengths must be positive"))
    return nothing
end

"""
    Detector(name, vertex, xarm, yarm, arm_lengths; check=true)
    Detector(name; latitude, longitude, elevation=0, xarm_azimuth,
             yarm_azimuth, xarm_tilt=0, yarm_tilt=0, arm_lengths, check=true)
    Detector(name)  # H1, L1, V1, K1

Immutable Earth-fixed interferometer geometry stored in `SVector`s. Cartesian
vectors have three components;
`xarm` and `yarm` are unit directions, and `arm_lengths` contains two full lengths in
metres. Nonorthogonal arms are supported. Geodetic coordinates use WGS-84, radians,
and metres; azimuth is clockwise from north and tilt positive above the horizontal.
Numeric inputs are promoted together without forcing Float64. Published built-in
directions retain their precision without normalization.

`check=false` bypasses numerical validation for AD/compiled calls with prevalidated
inputs. Vector sizes and numeric element types are always checked. Geometry is
physically real; the `Number` bound also permits Reactant's traced scalar types.
"""
struct Detector{T <: Number}
    name::String
    vertex::SVector{3, T}
    xarm::SVector{3, T}
    yarm::SVector{3, T}
    arm_lengths::SVector{2, T}

    function Detector{T}(
            name::String, vertex::SVector{3, T}, xarm::SVector{3, T},
            yarm::SVector{3, T}, arm_lengths::SVector{2, T}; check = true
        ) where {T <: Number}
        check && _check_geometry(vertex, xarm, yarm, arm_lengths)
        return new{T}(name, vertex, xarm, yarm, arm_lengths)
    end
end

function Detector(name::AbstractString, vertex, xarm, yarm, arm_lengths; check = true)
    length(vertex) == length(xarm) == length(yarm) == 3 ||
        throw(DimensionMismatch("Cartesian vectors must have three components"))
    length(arm_lengths) == 2 || throw(DimensionMismatch("arm_lengths must have two components"))
    p = promote(vertex..., xarm..., yarm..., arm_lengths...)
    return Detector{typeof(p[1])}(String(name), SVector(p[1:3]), SVector(p[4:6]), SVector(p[7:9]), SVector(p[10:11]); check)
end

function Detector(
        name::AbstractString; latitude = nothing, longitude = nothing, xarm_azimuth = nothing,
        yarm_azimuth = nothing, arm_lengths = nothing, elevation = nothing,
        xarm_tilt = nothing, yarm_tilt = nothing, check = true
    )
    coordinates = (latitude, longitude, xarm_azimuth, yarm_azimuth, arm_lengths)
    if all(isnothing, (coordinates..., elevation, xarm_tilt, yarm_tilt))
        key = Symbol(name)
        hasproperty(BUILTIN_GEOMETRY, key) || throw(ArgumentError("unknown detector: $name (choose H1, L1, V1, K1)"))
        return Detector(name, getproperty(BUILTIN_GEOMETRY, key)...; check)
    end
    any(isnothing, coordinates) && throw(ArgumentError("supply latitude, longitude, both azimuths, and arm_lengths"))
    elevation = something(elevation, 0)
    xarm_tilt = something(xarm_tilt, 0)
    yarm_tilt = something(yarm_tilt, 0)
    if check
        all(_finite_real, (latitude, longitude, elevation, xarm_azimuth, yarm_azimuth, xarm_tilt, yarm_tilt)) ||
            throw(ArgumentError("geodetic coordinates must be finite and real"))
        abs(latitude) <= π / 2 || throw(ArgumentError("latitude must be in [-π/2, π/2]"))
        all(t -> abs(t) <= π / 2, (xarm_tilt, yarm_tilt)) ||
            throw(ArgumentError("arm tilts must be in [-π/2, π/2]"))
    end
    length(arm_lengths) == 2 || throw(DimensionMismatch("arm_lengths must have two components"))
    latitude, longitude, elevation, xarm_azimuth, yarm_azimuth, xarm_tilt, yarm_tilt, lx, ly =
        promote(latitude, longitude, elevation, xarm_azimuth, yarm_azimuth, xarm_tilt, yarm_tilt, arm_lengths...)
    arm_lengths = (lx, ly)
    slat, clat = sincos(latitude)
    slon, clon = sincos(longitude)
    # WGS-84 equatorial radius and inverse flattening, in the input arithmetic.
    a = oftype(clat, 6378137)
    flattening = inv(oftype(clat, 298.257223563))
    e2 = flattening * (2 - flattening)
    radius = a / sqrt(1 - e2 * slat^2)
    vertex = (
        (radius + elevation) * clat * clon,
        (radius + elevation) * clat * slon,
        ((1 - e2) * radius + elevation) * slat,
    )
    north = (-slat * clon, -slat * slon, clat)
    east = (-slon, clon, zero(clat))
    up = (clat * clon, clat * slon, slat)
    arm(az, tilt) = ntuple(i -> cos(tilt) * (cos(az) * north[i] + sin(az) * east[i]) + sin(tilt) * up[i], 3)
    return Detector(name, vertex, arm(xarm_azimuth, xarm_tilt), arm(yarm_azimuth, yarm_tilt), arm_lengths; check)
end

"""
    detector_tensor(detector)

Return the symmetric response `(xarm ⊗ xarm - yarm ⊗ yarm)/2` as a 3×3 `SMatrix`.
Arm lengths do not enter the long-wavelength response.
"""
function detector_tensor(d::Detector)
    return (d.xarm * d.xarm' - d.yarm * d.yarm') / 2
end
