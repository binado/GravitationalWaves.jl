module GWDetectors

import GWSeries
using GWSeries: FrequencySeries

using LinearAlgebra: cross, dot
using StaticArrays: SVector, @SMatrix

export Detector, detectors, detector_tensor, greenwich_mean_sidereal_time,
    antenna_pattern, time_delay_from_geocenter, time_delay, light_travel_time, project_fd

const C_SI = 299792458

include("geometry.jl")
include("builtins.jl")
include("sidereal.jl")
include("response.jl")

end
