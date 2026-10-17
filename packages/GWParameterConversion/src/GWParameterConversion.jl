"""
    GWParameterConversion

Conversions between gravitational-wave source parameters (masses, spins, frames),
following bilby's names and conventions. Everything is generic over `Number`, so the
functions work with automatic differentiation, and the package has no dependencies.
Physical inputs must be real-valued; accepting Number accommodates traced real scalars.

Masses follow bilby's convention `mass_1 ≥ mass_2` and `mass_ratio = mass_2 / mass_1 ≤ 1`.
"""
module GWParameterConversion

export chirp_mass, total_mass, mass_ratio, symmetric_mass_ratio, reduced_mass
export total_mass_and_mass_ratio_to_component_masses,
    chirp_mass_and_mass_ratio_to_component_masses,
    chirp_mass_and_mass_ratio_to_total_mass,
    chirp_mass_and_total_mass_to_symmetric_mass_ratio,
    symmetric_mass_ratio_to_mass_ratio,
    mass_1_and_chirp_mass_to_mass_ratio,
    mass_2_and_chirp_mass_to_mass_ratio,
    source_to_detector_mass, detector_to_source_mass
export component_masses_to_chirp_mass_and_mass_ratio,
    component_masses_to_total_mass_and_mass_ratio
export effective_spin, effective_precessing_spin, is_aligned, aligned_spins_to_lal_frame,
    transform_precessing_initial_conditions, spins_to_lal_frame
export convert_parameters, InconsistentParametersError

include("constants.jl")
include("masses.jl")
include("mass_coordinates.jl")
include("spins.jl")
include("frames.jl")
include("convert.jl")

end
