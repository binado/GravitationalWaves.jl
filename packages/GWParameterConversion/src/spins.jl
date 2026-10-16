# Spin conversions ported from bilby/gw/conversion.py (`generate_spin_parameters`).

"""
    effective_spin(spin_1z, spin_2z, mass_ratio)

Mass-weighted aligned spin `χ_eff = (χ₁z + q χ₂z) / (1 + q)`.
"""
function effective_spin(spin_1z::Real, spin_2z::Real, mass_ratio::Real)
    return (spin_1z + spin_2z * mass_ratio) / (1 + mass_ratio)
end

"""
    effective_precessing_spin(spin_1x, spin_1y, spin_2x, spin_2y, mass_ratio)

Precession parameter `χ_p = max(χ₁⊥, (4q+3)/(3q+4) q χ₂⊥)` from the in-plane spin components.
"""
function effective_precessing_spin(
        spin_1x::Real, spin_1y::Real, spin_2x::Real, spin_2y::Real, mass_ratio::Real
    )
    chi_1_in_plane = sqrt(spin_1x^2 + spin_1y^2)
    chi_2_in_plane = sqrt(spin_2x^2 + spin_2y^2)
    return effective_precessing_spin(chi_1_in_plane, chi_2_in_plane, mass_ratio)
end

"""
    effective_precessing_spin(chi_1_in_plane, chi_2_in_plane, mass_ratio)

Same from the in-plane spin magnitudes.
"""
function effective_precessing_spin(
        chi_1_in_plane::Real, chi_2_in_plane::Real, mass_ratio::Real
    )
    return max(
        chi_1_in_plane,
        (4 * mass_ratio + 3) / (3 * mass_ratio + 4) * mass_ratio * chi_2_in_plane,
    )
end
