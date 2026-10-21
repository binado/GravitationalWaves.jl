"""Read one-sided spectral densities and sample stationary Gaussian strain noise."""
module GWNoise

import GWTimeFrequencyDomains

using GWTimeFrequencyDomains: TimeDomain, FrequencyDomain, nsamples, sampling_frequency, channels
using GWTimeFrequencyDomains: epoch, frequency_spacing, check_coordinate, check_same_frequency_grid
using GWTimeFrequencyDomains: to_timedomain
import GWTimeFrequencyDomains: frequencies, check_compatible

using DelimitedFiles: readdlm
# Activate GWTimeFrequencyDomains' FFT extension.
import FFTW
import LinearAlgebra
using Random: AbstractRNG, default_rng, randn

export available_psds, psd_path
export TabulatedPSD, read_psd, read_asd, psd, asd
export GaussianNoise, frequencies, noise_fd, noise_fd!, noise_td

include("psd.jl")
include("bundled.jl")
include("noise.jl")

end
