"""Shared contract, reference and gradient test suites for `AbstractWaveforms` generators."""
module AbstractWaveformsTests

using Test
using AbstractWaveforms
using AbstractWaveforms: AlignedSpinQuadrupolarWaveformGenerator, amplitude_phase_fd, check_frequency_grid
using GWSeries
using GWMatchedFiltering
using GWMatchedFiltering: FrequencyDomainMetric, complex_inner_product, optimal_snr_squared
using ADTypes: AutoFiniteDifferences
using DifferentiationInterface: gradient, jacobian, prepare_gradient, value_and_gradient
using FiniteDifferences: central_fdm

export test_generator_contract, test_amplitude_phase_fd_contract, test_modes_td_contract, test_modes_fd_contract
export load_cases, default_cases, applicable_cases, fd_parameters, surrogate_parameters
export scaled_coordinates, fd_adapter, surrogate_adapter
export mismatch, flat_metric, test_reference
export FixedObjective, test_gradients, test_snr_gradients, test_phase_gradients, test_waveform_jacobian

include("contract.jl")
include("cases.jl")
include("mismatch.jl")
include("reference.jl")
include("gradients.jl")

end
