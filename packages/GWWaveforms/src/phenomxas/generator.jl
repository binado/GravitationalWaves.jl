"""
    IMRPhenomXAS(; reference_frequency=50.0, minimum_frequency=20.0,
                  maximum_frequency=nothing, pn_amplitude_order=0)

Default aligned-spin, quadrupole IMRPhenomXAS model (Pratten et al. 2020).
Frequencies are in Hz; `maximum_frequency=nothing` uses the grid endpoint.
The intrinsic cutoff is `0.3 / M_sec`. No tides, precession, higher modes, or
alternative LAL calibration-version flags are supported.
"""
struct IMRPhenomXAS <: AlignedSpinQuadrupolarWaveformGenerator
    reference_frequency::Float64
    minimum_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    pn_amplitude_order::Int
end

function IMRPhenomXAS(;
        reference_frequency::Real = 50.0, minimum_frequency::Real = 20.0,
        maximum_frequency::Union{Nothing, Real} = nothing, pn_amplitude_order::Integer = 0,
    )
    check_pn_orders(; pn_amplitude_order)
    return IMRPhenomXAS(
        Float64(reference_frequency), Float64(minimum_frequency),
        maximum_or_nothing(maximum_frequency), Int(pn_amplitude_order)
    )
end

approximant(::IMRPhenomXAS) = "IMRPhenomXAS"
reference_frequency(g::IMRPhenomXAS) = g.reference_frequency
minimum_frequency(g::IMRPhenomXAS) = g.minimum_frequency
maximum_frequency(g::IMRPhenomXAS) = g.maximum_frequency
spin_support(::IMRPhenomXAS) = :aligned
is_differentiable(::IMRPhenomXAS) = true
