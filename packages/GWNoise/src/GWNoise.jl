"""Read one-sided spectral densities and sample stationary Gaussian strain noise."""
module GWNoise

import GWTimeFrequencyDomains

using GWTimeFrequencyDomains: TimeDomain, FrequencyDomain, nsamples, sampling_frequency, channels
using GWTimeFrequencyDomains: epoch, frequency_spacing, check_coordinate, check_same_frequency_grid
import GWTimeFrequencyDomains: frequencies, check_compatible

using DelimitedFiles: readdlm
using FFTW: irfft
using Random: AbstractRNG, default_rng, randn

export TabulatedPSD, read_psd, read_asd, psd, asd
export GaussianNoise, frequencies, noise_fd, noise_fd!, noise_td, noise_td!

include("psd.jl")
include("noise.jl")

end
