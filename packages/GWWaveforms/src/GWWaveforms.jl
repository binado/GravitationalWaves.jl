"""
    GWWaveforms

Pure-Julia frequency-domain compact-binary waveforms, generic over `Real` and so
differentiable with ForwardDiff, Enzyme and Mooncake. Each model is a generator type
(`TaylorF2`, `IMRPhenomD`) implementing the `AbstractGWWaveforms` interface (`hphc_fd`, bilby's
parameter names) shared with `LALWaveforms`. The implementations are ports of LALSimulation's
TaylorF2 and IMRPhenomD, which `LALWaveforms` reproduces for testing.
"""
module GWWaveforms

using GWParameterConversion: MTSUN_SI, PC_SI, MASS_KEYS, convert_parameters,
    is_aligned, aligned_spins_to_lal_frame

using AbstractGWWaveforms: AbstractWaveformGenerator, check_frequency_grid, frequency_bounds
import AbstractGWWaveforms: hphc_fd, approximant, reference_frequency, minimum_frequency,
    maximum_frequency, spin_support, supports_mode_array, is_differentiable

export AbstractWaveformGenerator, hphc_fd
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable
export TaylorF2, IMRPhenomD

# Building blocks for later interfaces: one-frequency evaluation of the models.
public taylorf2_pn_phasing, taylorf2_phasing, phenomd_model, phenomd_amplitude, phenomd_phase

include("constants.jl")
include("generators.jl")
include("parameters.jl")
include("taylorf2.jl")
include("phenomd.jl")
include("fd.jl")

end
