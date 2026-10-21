"""Read one-sided spectral densities and sample stationary Gaussian strain noise."""
module GWNoise

import GWSeries

using GWSeries: TimeSeries, FrequencySeries, metadata, sampling_interval, nsamples, sampling_frequency, channels
using GWSeries: epoch, frequency_spacing, check_coordinate, check_same_frequency_grid
using GWSeries: to_timedomain
import GWSeries: frequencies, check_compatible

using DelimitedFiles: readdlm
using Random: AbstractRNG, default_rng, randn

export available_psds, psd_path
export TabulatedPSD, read_psd, read_asd, psd, asd
export GaussianNoise, frequencies, noise_fd, noise_fd!, noise_td

include("psd.jl")
include("bundled.jl")
include("noise.jl")

end
