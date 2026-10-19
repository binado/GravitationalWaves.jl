"""
    Approximant

Supertype of the waveform models implemented here, each a concrete
[`AbstractWaveformGenerator`](@ref) carrying only the options that model uses. Defaults are
those of bilby's `lal_binary_black_hole`.
"""
abstract type Approximant <: AbstractWaveformGenerator end

function check_pn_orders(; pn_spin_order = 0, pn_phase_order = 0, pn_amplitude_order = 0)
    -1 <= pn_spin_order <= 7 || throw(ArgumentError("invalid pn_spin_order $pn_spin_order"))
    -1 <= pn_phase_order <= 7 || throw(ArgumentError("invalid pn_phase_order $pn_phase_order"))
    # LAL has no 1PN-relative amplitude correction.
    (pn_amplitude_order == -1 || pn_amplitude_order == 0 || 2 <= pn_amplitude_order <= 7) ||
        throw(ArgumentError("invalid pn_amplitude_order $pn_amplitude_order"))
    return nothing
end

maximum_or_nothing(f::Nothing) = nothing
maximum_or_nothing(f::Real) = Float64(f)

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
struct TaylorF2 <: Approximant
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

"""
    IMRPhenomD(; kwargs...)

Aligned-spin inspiral-merger-ringdown model of Husa et al. and Khan et al.
(`XLALSimIMRPhenomDGenerateFD`).

# Keywords
- `reference_frequency = 50.0`, `minimum_frequency = 20.0` (Hz)
- `maximum_frequency = nothing`: the last frequency of the grid
- `pn_amplitude_order = 0`: only shifts the start frequency, as in LAL.
"""
struct IMRPhenomD <: Approximant
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

"""
    IMRPhenomXAS(; reference_frequency=50.0, minimum_frequency=20.0,
                  maximum_frequency=nothing, pn_amplitude_order=0)

Default aligned-spin, quadrupole IMRPhenomXAS model (Pratten et al. 2020).
Frequencies are in Hz; `maximum_frequency=nothing` uses the grid endpoint.
The intrinsic cutoff is `0.3 / M_sec`. No tides, precession, higher modes, or
alternative LAL calibration-version flags are supported.
"""
struct IMRPhenomXAS <: Approximant
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

approximant(::A) where {A <: Approximant} = String(nameof(A))
reference_frequency(g::Approximant) = g.reference_frequency
minimum_frequency(g::Approximant) = g.minimum_frequency
maximum_frequency(g::Approximant) = g.maximum_frequency

spin_support(::Approximant) = :aligned
is_differentiable(::Approximant) = true
