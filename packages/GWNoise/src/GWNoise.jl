"""Read one-sided spectral densities and sample stationary Gaussian strain noise."""
module GWNoise

import GWSeries

using GWSeries: TimeSeries, FrequencySeries, metadata, sampling_interval, nsamples, channels
using GWSeries: check_coordinate, check_same_frequency_grid, time_samples, to_timedomain
import GWSeries: frequencies, frequency_spacing, check_compatible

using DelimitedFiles: readdlm
using LinearAlgebra: Hermitian, I, PosDefException, cholesky, ishermitian
using Random: AbstractRNG, default_rng, randn

export available_psds, psd_path
export TabulatedPSD, read_psd, read_asd, psd, asd
export band_psd, whitening_factor, sample_noise, sample_noise!
export GaussianNoise, frequencies, noise_fd, noise_fd!, noise_td

include("psd.jl")
include("bundled.jl")
include("noise.jl")
include("arrays.jl")

end
