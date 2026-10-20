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

