"""
    GWWaveforms

Pure-Julia frequency-domain compact-binary waveforms, generic over `Real` and so
differentiable with ForwardDiff and Mooncake. Each model is a generator type
(`TaylorF2`, `IMRPhenomD`, `IMRPhenomXAS`) implementing the `AbstractGWWaveforms` interface (`hphc_fd`, bilby's
parameter names) shared with `LALWaveforms`. The implementations are ports of LALSimulation's
TaylorF2, IMRPhenomD and IMRPhenomXAS, with `LALWaveforms` as the testing oracle.
"""
module GWWaveforms

using StaticArrays: SMatrix, SVector
using LinearAlgebra: lu

using GWParameterConversion: MTSUN_SI, PC_SI, MASS_KEYS, convert_parameters,
    is_aligned, aligned_spins_to_lal_frame

using AbstractGWWaveforms: AbstractWaveformGenerator, AlignedSpinQuadrupolarWaveformGenerator,
    check_frequency_grid, frequency_bounds
import AbstractGWWaveforms: hphc_fd, amplitude_phase_fd, validate_parameters, approximant, reference_frequency, minimum_frequency,
    maximum_frequency, spin_support, supports_mode_array, is_differentiable

export AbstractWaveformGenerator, AlignedSpinQuadrupolarWaveformGenerator
export hphc_fd, amplitude_phase_fd, validate_parameters
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable
export TaylorF2, IMRPhenomD, IMRPhenomXAS

# Building blocks for later interfaces: one-frequency evaluation of the models.
public taylorf2_pn_phasing, taylorf2_phasing, phenomd_model, phenomd_amplitude, phenomd_phase

include("constants.jl")
include("generators.jl")
include("parameters.jl")
include("taylorf2.jl")
include("phenom_common.jl")
include("phenomd/model.jl")
include("phenomxas/fits.jl")
include("phenomxas/utils.jl")
include("phenomxas/model.jl")
include("fd.jl")

end
