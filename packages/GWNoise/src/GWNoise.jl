"""Read one-sided spectral densities and sample stationary Gaussian noise on plain arrays."""
module GWNoise

using DelimitedFiles: readdlm
using LinearAlgebra: Hermitian, I, PosDefException, cholesky, ishermitian
using Random: AbstractRNG, default_rng, randn

export available_psds, psd_path
export TabulatedPSD, psd, asd
export band_psd, whitening_factor, sample_noise, sample_noise!

include("psd.jl")
include("bundled.jl")
include("arrays.jl")

end
