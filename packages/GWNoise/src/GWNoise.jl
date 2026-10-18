"""Read one-sided spectral densities and sample stationary Gaussian strain noise."""
module GWNoise

using GWSeries: TimeSeries, FrequencySeries, epoch, channels, check_same_grid
import GWSeries: frequencies

using DelimitedFiles: readdlm
using FFTW: irfft
using Random: AbstractRNG, default_rng, randn

export TabulatedPSD, read_psd, read_asd, psd, asd
export GaussianNoise, frequencies, noise_fd, noise_fd!, noise_td, noise_td!

include("psd.jl")
include("noise.jl")

end
