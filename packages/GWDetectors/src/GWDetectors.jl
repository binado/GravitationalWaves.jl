module GWDetectors

import TOML
import GWSeries
using GWSeries: FrequencySeries, metadata, channels, frequencies, frequency_spacing, nsamples, epoch

using LinearAlgebra: cross, dot
using StaticArrays: SVector, @SMatrix

export Detector, detectors, sensitivities, detector_tensor, greenwich_mean_sidereal_time,
    antenna_pattern, time_delay_from_geocenter, time_delay, light_travel_time, project_fd

const C_SI = 299792458

include("geometry.jl")
include("builtins.jl")
include("sensitivity.jl")
include("sidereal.jl")
include("response.jl")

end
