"""
    LALWaveforms

Frequency-domain compact-binary waveforms from LALSimulation, called directly
through `ccall` with bilby's parameter names, units and parameter conversions.

The LAL shared libraries are not bundled: point the package at them with the
`LALWAVEFORMS_LIBDIR` environment variable, [`set_libdir!`](@ref) or
[`set_libdir_from_python!`](@ref).
"""
module LALWaveforms

using Libdl
using Preferences

export hphc_fd, WaveformSettings, LALError, set_libdir!, set_libdir_from_python!

include("types.jl")
include("library.jl")
include("errors.jl")
include("settings.jl")
include("parameters.jl")
include("conversion.jl")
include("fd.jl")

end
