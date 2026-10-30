"""Uniform sampled DimensionalData matrices with reconstructible grids and physical AbstractFFTs transforms."""
module GWSeries

import AbstractFFTs
import DimensionalData as DD
using DimensionalData: AbstractDimArray, DimArray, Dimension, TimeDim, @dim, dims, dimnum, lookup, rebuild
using DimensionalData: Dimensions, NoName, NoMetadata, At, Near, var".."
using DimensionalData.Lookups: Sampled, Categorical, ForwardOrdered, Unordered, Regular, Points
using LinearAlgebra: mul!, Transpose, Adjoint

export AbstractSeries, TimeSeries, FrequencySeries, grid_metadata, flatten
export GPSTime, Freq, Chan, TimeGrid, SpectralGrid
# DimensionalData selectors, so `x[Chan = At(:H1)]` needs no extra import.
export At, Near, var".."
export epoch, channels, nsamples, times, frequencies, sampling_interval, sampling_frequency, duration, frequency_spacing, time_samples
export check_data, check_same_grid, check_same_frequency_grid, check_compatible, check_fourier_domain
export frequency_bins, crop, pad, to_frequencydomain, to_timedomain
export plan_to_frequencydomain, plan_to_timedomain, transform, transform!
export PlanckWindow, apply_window, apply_window!

# Construction helpers that dependent packages build prepared objects from.
public time_metadata, frequency_metadata, metadata_value, channel_order, check_coordinate, ChannelNames

include("grid_common.jl")
include("grids.jl")
include("dims.jl")
include("series.jl")
include("grid_operations.jl")
include("windows.jl")
include("fourier.jl")

end
