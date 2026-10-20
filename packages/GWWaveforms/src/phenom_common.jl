# One-based indices into the requested grid; no relationship between index and frequency
# is assumed. Grid and generator validation runs separately on the host.
function phenom_frequency_slice(f::AbstractRange{Float64}, f_start, f_min, f_max)
    return (searchsortedfirst(f, max(f_start, f_min)), searchsortedlast(f, f_max))
end

# LAL uses the starting frequency when no reference frequency is specified.
phenom_reference_frequency(g, f_start) = g.reference_frequency == 0.0 ? f_start : g.reference_frequency

# Shared strain scale; keep LAL's multiplication order for numerical fidelity.
function phenom_amplitude_scale(m1, m2, luminosity_distance)
    M = m1 + m2
    distance = luminosity_distance * 1.0e6 * PC_SI
    return 2.0 * sqrt(5.0 / (64.0 * Float64(π))) * M * MRSUN_SI * M * MTSUN_SI / distance
end

# Host-only physical checks; remnant and frequency restrictions belong to each model.
function check_phenom_parameters(params::NamedTuple)
    (params.mass_1 > 0 && params.mass_2 > 0) || throw(ArgumentError("masses must be positive"))
    params.luminosity_distance > 0 || throw(ArgumentError("luminosity distance must be positive"))
    (abs(params.chi_1) <= 1 && abs(params.chi_2) <= 1) ||
        throw(ArgumentError("spins outside the range [-1, 1] are not supported"))
    return nothing
end
