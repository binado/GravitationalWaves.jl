"""
    TaylorF2(; kwargs...)

Post-Newtonian stationary-phase inspiral (`XLALSimInspiralTaylorF2`), aligned spins only.

# Keywords
- `reference_frequency = 50.0`, `minimum_frequency = 20.0` (Hz)
- `maximum_frequency = nothing`: the last frequency of the grid
- `pn_spin_order = -1`, `pn_phase_order = -1`, `pn_amplitude_order = 0`; `-1` selects the highest
  implemented order. The amplitude order is the SPA amplitude correction and shifts the start
  frequency.

Tidal parameters are not supported.
"""
struct TaylorF2 <: AlignedSpinQuadrupolarWaveformGenerator
    reference_frequency::Float64
    minimum_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    pn_spin_order::Int
    pn_phase_order::Int
    pn_amplitude_order::Int
end

function TaylorF2(;
        reference_frequency::Real = 50.0,
        minimum_frequency::Real = 20.0,
        maximum_frequency::Union{Nothing, Real} = nothing,
        pn_spin_order::Integer = -1,
        pn_phase_order::Integer = -1,
        pn_amplitude_order::Integer = 0,
    )
    check_pn_orders(; pn_spin_order, pn_phase_order, pn_amplitude_order)
    return TaylorF2(
        Float64(reference_frequency), Float64(minimum_frequency),
        maximum_or_nothing(maximum_frequency),
        Int(pn_spin_order), Int(pn_phase_order), Int(pn_amplitude_order),
    )
end

approximant(::TaylorF2) = "TaylorF2"
reference_frequency(g::TaylorF2) = g.reference_frequency
minimum_frequency(g::TaylorF2) = g.minimum_frequency
maximum_frequency(g::TaylorF2) = g.maximum_frequency
spin_support(::TaylorF2) = :aligned
is_differentiable(::TaylorF2) = true
