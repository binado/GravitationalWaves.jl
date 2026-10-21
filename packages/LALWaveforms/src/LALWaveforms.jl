"""
    LALWaveforms

Frequency-domain compact-binary waveforms from LALSimulation, called directly
through `ccall` with bilby's parameter names, units and parameter conversions. A
[`LALWaveformGenerator`](@ref) implements the `AbstractWaveforms` interface shared with
`FrequencyDomainWaveforms`.

The LAL shared libraries are not bundled: point the package at them with the
`LALWAVEFORMS_LIBDIR` environment variable, [`set_libdir!`](@ref) or
[`set_libdir_from_python!`](@ref).
"""
module LALWaveforms

using GWParameterConversion: MSUN_SI, PC_SI, MASS_KEYS, convert_parameters, is_aligned,
    aligned_spins_to_lal_frame
using AbstractWaveforms: AbstractWaveformGenerator, check_frequency_grid, frequency_bounds
import AbstractWaveforms: hphc_fd, validate_parameters, approximant, reference_frequency, minimum_frequency,
    maximum_frequency, spin_support, supports_mode_array, is_differentiable
using Libdl
using Preferences

export LALWaveformGenerator, LALError, set_libdir!, set_libdir_from_python!
export AbstractWaveformGenerator, hphc_fd, validate_parameters
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable

include("types.jl")
include("library.jl")
include("errors.jl")
include("generator.jl")
include("parameters.jl")
include("conversion.jl")
include("fd.jl")

end
