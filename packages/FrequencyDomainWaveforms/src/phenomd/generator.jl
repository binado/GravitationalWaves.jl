"""
    IMRPhenomD(; kwargs...)

Aligned-spin inspiral-merger-ringdown model of Husa et al. and Khan et al.
(`XLALSimIMRPhenomDGenerateFD`).

# Keywords
- `reference_frequency = 50.0`, `minimum_frequency = 20.0` (Hz)
- `maximum_frequency = nothing`: the last frequency of the grid
- `pn_amplitude_order = 0`: only shifts the start frequency, as in LAL.
"""
struct IMRPhenomD <: AlignedSpinQuadrupolarWaveformGenerator
    reference_frequency::Float64
    minimum_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    pn_amplitude_order::Int
end

function IMRPhenomD(;
        reference_frequency::Real = 50.0,
        minimum_frequency::Real = 20.0,
        maximum_frequency::Union{Nothing, Real} = nothing,
        pn_amplitude_order::Integer = 0,
    )
    check_pn_orders(; pn_amplitude_order)
    return IMRPhenomD(
        Float64(reference_frequency), Float64(minimum_frequency),
        maximum_or_nothing(maximum_frequency), Int(pn_amplitude_order),
    )
end

approximant(::IMRPhenomD) = "IMRPhenomD"
reference_frequency(g::IMRPhenomD) = g.reference_frequency
minimum_frequency(g::IMRPhenomD) = g.minimum_frequency
maximum_frequency(g::IMRPhenomD) = g.maximum_frequency
spin_support(::IMRPhenomD) = :aligned
is_differentiable(::IMRPhenomD) = true
