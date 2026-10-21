"""Uniform sampled matrices with reconstructible metadata and physical AbstractFFTs transforms."""
module GWSeries

import AbstractFFTs
using LinearAlgebra: mul!, Transpose, Adjoint

export AbstractSeries, TimeSeries, FrequencySeries, metadata, flatten
export epoch, channels, nsamples, times, frequencies, sampling_interval, sampling_frequency, duration, frequency_spacing, time_samples
export check_data, check_same_grid, check_same_frequency_grid, check_compatible, check_fourier_domain
export frequency_bins, crop, pad, to_frequencydomain, to_timedomain
export plan_to_frequencydomain, plan_to_timedomain, transform, transform!
export PlanckWindow, apply_window, apply_window!

# Construction helpers that dependent packages build prepared objects from.
public time_metadata, frequency_metadata, metadata_value, channel_order, check_coordinate

include("grid_common.jl")
include("series.jl")
include("indexing.jl")
include("arithmetic.jl")
include("grid_operations.jl")
include("windows.jl")
include("fourier.jl")

end
