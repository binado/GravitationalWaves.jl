"""
    GWWaveforms

Pure-Julia frequency-domain compact-binary waveforms, generic over `Real` and so
differentiable with ForwardDiff, Enzyme and Mooncake. The API mirrors `LALWaveforms`
(`WaveformSettings`, `hphc_fd`, bilby's parameter names), and the implementations are ports of
LALSimulation's TaylorF2 and IMRPhenomD, which `LALWaveforms` reproduces for testing.
"""
module GWWaveforms

using GWParameterConversion: MSUN_SI, MTSUN_SI, PC_SI, MASS_KEYS, convert_parameters,
    is_aligned, aligned_spins_to_lal_frame

export hphc_fd, WaveformSettings

# Building blocks for later interfaces: one-frequency evaluation of the models.
public taylorf2_pn_phasing, taylorf2_phasing, phenomd_model, phenomd_amplitude, phenomd_phase

include("constants.jl")
include("approximants.jl")
include("settings.jl")
include("parameters.jl")
include("taylorf2.jl")
include("phenomd.jl")
include("fd.jl")

end
