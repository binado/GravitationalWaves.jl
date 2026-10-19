module GWSEOBNR

using AbstractGWWaveforms
using GWParameterConversion
using ForwardDiff
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using FFTW
using OrdinaryDiffEqVerner: Vern9
using SciMLBase

import AbstractGWWaveforms: hphc_fd, validate_parameters, approximant,
    reference_frequency, minimum_frequency, maximum_frequency,
    spin_support, supports_mode_array, is_differentiable
using AbstractGWWaveforms: check_frequency_grid, frequency_bounds

export SEOBNRv5HM, hphc_fd, validate_parameters
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable

include("generator.jl")
include("parameters.jl")
include("fits.jl")
include("hamiltonian.jl")
include("fft.jl")
include("tail.jl")
include("flux_coefficients.jl")
include("flux.jl")
include("mode_fits.jl")
include("remnant.jl")
include("roots.jl")
include("dynamics.jl")
include("qnm_data.jl")
include("qnm.jl")
include("merger.jl")
include("fd.jl")

end
