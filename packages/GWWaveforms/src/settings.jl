"""
    WaveformSettings(approximant; kwargs...)

Everything about a waveform call that does not vary between parameter draws: the
approximant, frequency bounds and post-Newtonian orders. Mirrors
`LALWaveforms.WaveformSettings`, so code can switch between the LAL-backed and the native
implementation. Defaults are those of bilby's `lal_binary_black_hole`.

`approximant` is `"TaylorF2"` or `"IMRPhenomD"` (a `Symbol` also works); other names throw
`ArgumentError`.

# Keywords
- `reference_frequency = 50.0`, `minimum_frequency = 20.0` (Hz)
- `maximum_frequency = nothing`: the last frequency of the grid
- `pn_spin_order = -1`, `pn_tidal_order = -1`, `pn_phase_order = -1`,
  `pn_amplitude_order = 0`; `-1` selects the highest implemented order. The spin and phase
  orders only act on TaylorF2, the amplitude order is TaylorF2's SPA amplitude correction and
  shifts the start frequency for both models. Tidal parameters are not supported, so the tidal
  order has no effect.
- `mode_array = nothing`: neither model has higher modes, so anything else throws
  `ArgumentError`.
"""
struct WaveformSettings{A <: Approximant}
    approximant::A
    reference_frequency::Float64
    minimum_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    pn_spin_order::Int
    pn_tidal_order::Int
    pn_phase_order::Int
    pn_amplitude_order::Int
end

function WaveformSettings(
        approximant::Union{Symbol, AbstractString};
        reference_frequency::Real = 50.0,
        minimum_frequency::Real = 20.0,
        maximum_frequency::Union{Nothing, Real} = nothing,
        pn_spin_order::Integer = -1,
        pn_tidal_order::Integer = -1,
        pn_phase_order::Integer = -1,
        pn_amplitude_order::Integer = 0,
        mode_array = nothing,
    )
    model = approximant_from_name(approximant)
    mode_array === nothing || throw(
        ArgumentError("`mode_array` is not supported: $(approximant_name(model)) has no higher modes")
    )
    -1 <= pn_spin_order <= 7 || throw(ArgumentError("invalid pn_spin_order $pn_spin_order"))
    -1 <= pn_phase_order <= 7 || throw(ArgumentError("invalid pn_phase_order $pn_phase_order"))
    # LAL has no 1PN-relative amplitude correction.
    (pn_amplitude_order == -1 || pn_amplitude_order == 0 || 2 <= pn_amplitude_order <= 7) ||
        throw(ArgumentError("invalid pn_amplitude_order $pn_amplitude_order"))
    return WaveformSettings(
        model, Float64(reference_frequency), Float64(minimum_frequency),
        maximum_frequency === nothing ? nothing : Float64(maximum_frequency),
        Int(pn_spin_order), Int(pn_tidal_order), Int(pn_phase_order), Int(pn_amplitude_order),
    )
end

function Base.show(io::IO, s::WaveformSettings)
    return print(
        io, "WaveformSettings(", repr(approximant_name(s.approximant)),
        "; reference_frequency=", s.reference_frequency,
        ", minimum_frequency=", s.minimum_frequency,
        ", maximum_frequency=", s.maximum_frequency, ")"
    )
end
