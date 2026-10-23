const FIELDS = (
    :mass_1_msun, :mass_2_msun, :spin_1x, :spin_1y, :spin_1z,
    :spin_2x, :spin_2y, :spin_2z, :distance_mpc, :inclination_rad, :phase_rad,
)

"""
    load_cases(path = <bundled parameters.csv>)

Read the shared physical-parameter table as `Dict(id => NamedTuple)`. Masses are solar masses,
distances Mpc, angles radians, and spins dimensionless Cartesian components in the LAL source frame.
"""
function load_cases(path = joinpath(@__DIR__, "..", "data", "parameters.csv"))
    lines = readlines(path)
    split(first(lines), ',') == ["id"; string.(FIELDS)...] || error("unexpected waveform CSV schema")
    cases = Dict{String, NamedTuple{FIELDS, NTuple{11, Float64}}}()
    for line in lines[2:end]
        row = split(line, ',')
        length(row) == 12 || error("expected 12 columns")
        id = String(row[1])
        isempty(id) && error("empty case ID")
        haskey(cases, id) && error("duplicate case ID: $id")
        p = NamedTuple{FIELDS}(Tuple(parse.(Float64, row[2:end])))
        all(isfinite, p) || error("nonfinite case: $id")
        p.mass_1_msun >= p.mass_2_msun > 0 || error("unordered masses: $id")
        p.distance_mpc > 0 || error("nonpositive distance: $id")
        sum(abs2, (p.spin_1x, p.spin_1y, p.spin_1z)) <= 1 || error("invalid spin: $id")
        sum(abs2, (p.spin_2x, p.spin_2y, p.spin_2z)) <= 1 || error("invalid spin: $id")
        cases[id] = p
    end
    return cases
end

const DEFAULT_CASES = Ref{Dict{String, NamedTuple{FIELDS, NTuple{11, Float64}}}}()

"""The bundled cases, loaded once on first use."""
function default_cases()
    isassigned(DEFAULT_CASES) || (DEFAULT_CASES[] = load_cases())
    return DEFAULT_CASES[]
end

has_transverse_spin(p) = !all(iszero, (p.spin_1x, p.spin_1y, p.spin_2x, p.spin_2y))
has_spin(p) = has_transverse_spin(p) || !iszero(p.spin_1z) || !iszero(p.spin_2z)

function supports_case(g, p)
    support = spin_support(g)
    support === :precessing && return true
    support === :aligned && return !has_transverse_spin(p)
    return !has_spin(p)
end

"""
    applicable_cases(generators..., cases = default_cases())

Sorted `(id, case)` pairs that every given generator supports, judged by `spin_support`:
`:none` keeps only non-spinning rows, `:aligned` drops rows with transverse spin.
"""
applicable_cases(g::AbstractWaveformGenerator, gs...) = applicable_cases((g, gs...), default_cases())
applicable_cases(g::AbstractWaveformGenerator, cases::AbstractDict) = applicable_cases((g,), cases)
function applicable_cases(gs::Tuple, cases::AbstractDict = default_cases())
    return [(id, p) for (id, p) in sort!(collect(cases); by = first) if all(g -> supports_case(g, p), gs)]
end

"""Aligned-spin parameters with bilby names; throws for rows with transverse spin."""
function fd_parameters(p)
    has_transverse_spin(p) && throw(ArgumentError("FD fixtures must have no transverse spins"))
    return (
        mass_1 = p.mass_1_msun, mass_2 = p.mass_2_msun,
        a_1 = abs(p.spin_1z), a_2 = abs(p.spin_2z),
        tilt_1 = p.spin_1z < 0 ? Float64(π) : 0.0,
        tilt_2 = p.spin_2z < 0 ? Float64(π) : 0.0,
        luminosity_distance = p.distance_mpc, theta_jn = p.inclination_rad, phase = p.phase_rad,
    )
end

"""Surrogate parameters: total mass, q = m1/m2, Cartesian spins, inclination and reference phase."""
function surrogate_parameters(p)
    return (
        total_mass = p.mass_1_msun + p.mass_2_msun,
        mass_ratio = p.mass_1_msun / p.mass_2_msun,
        spin_1 = (p.spin_1x, p.spin_1y, p.spin_1z),
        spin_2 = (p.spin_2x, p.spin_2y, p.spin_2z),
        luminosity_distance = p.distance_mpc, inclination = p.inclination_rad, phi_ref = p.phase_rad,
    )
end

"""
    scaled_coordinates(p; precessing = false) -> (x, scales)

Dimensionless coordinates for differentiation: masses and distance are divided by their fixture
values, spins and angles are untouched. This avoids polar spin-coordinate singularities at zero spin.
"""
function scaled_coordinates(p; precessing = false)
    x = precessing ? collect(values(p)) :
        [p.mass_1_msun, p.mass_2_msun, p.spin_1z, p.spin_2z, p.distance_mpc, p.inclination_rad, p.phase_rad]
    scales = ones(length(x))
    scales[1:2] .= x[1:2]
    scales[end - 2] = p.distance_mpc
    return x ./ scales, scales
end

"""Map scaled aligned-spin coordinates to bilby parameters. Tilt is fixed at preparation so the differentiated amplitude is smooth through zero."""
function fd_adapter(p, scales)
    s1 = p.spin_1z < 0 ? -1.0 : 1.0
    s2 = p.spin_2z < 0 ? -1.0 : 1.0
    return z -> begin
        x = z .* scales
        (
            mass_1 = x[1], mass_2 = x[2], a_1 = s1 * x[3], a_2 = s2 * x[4],
            tilt_1 = s1 < 0 ? Float64(π) : 0.0, tilt_2 = s2 < 0 ? Float64(π) : 0.0,
            luminosity_distance = x[5], theta_jn = x[6], phase = x[7],
        )
    end
end
"""Map scaled full-spin coordinates to surrogate parameters."""
surrogate_adapter(scales) = z -> surrogate_parameters(NamedTuple{FIELDS}(Tuple(z .* scales)))
