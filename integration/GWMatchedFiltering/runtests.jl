using Test
using GWMatchedFiltering
using GWSeries
using GWDetectors
using GWWaveforms
using GWNoise
using LinearAlgebra
using ForwardDiff
using DifferentiationInterface
using Enzyme: Enzyme
using Mooncake: Mooncake
using Reactant

Reactant.set_default_backend("cpu")
include("../../packages/GWMatchedFiltering/test/ad_helpers.jl")

const PIPE_F = 0.0:1.0:128.0
# Projection accepts arbitrary vectors. Materialize its fixed grid once so AD
# does not enter Base's TwicePrecision range-broadcast machinery.
const PIPE_PROJECT_F = collect(PIPE_F)
const PIPE_GENERATOR = TaylorF2(; minimum_frequency = 20.0, maximum_frequency = 128.0)
const PIPE_SKY = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
const PIPE_METRIC = FrequencyDomainMetric(
    PIPE_F, (:H1, :L1),
    repeat(reshape([2.0 0.2; 0.2 1.0] .* 1.0e-46, 1, 2, 2), length(PIPE_F), 1, 1);
    minimum_frequency = 20.0,
)

function pipeline_signal(x)
    params = (; mass_1 = x[1], mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = x[2])
    pol = hphc_fd(PIPE_F, params, PIPE_GENERATOR)
    h1 = project_fd(PIPE_PROJECT_F, pol, PIPE_SKY, detectors.H1; start_time = 1126259460.0)
    l1 = project_fd(PIPE_PROJECT_F, pol, PIPE_SKY, detectors.L1; start_time = 1126259460.0)
    return hcat(h1, l1)
end

const PIPE_DATA = pipeline_signal([36.0, 1.3])
pipeline_loss(x) = log_likelihood_ratio(pipeline_signal(x), PIPE_DATA, PIPE_METRIC)

@testset "GWMatchedFiltering AD integration" begin
    for (loss, sample) in ((diagonal_loss, [0.8, 0.2]), (correlated_loss, [0.8, 0.2]), (statistics_loss, [0.8, 0.2]), (pipeline_loss, [35.8, 1.2]))
        @info "Checking GWMatchedFiltering derivatives" objective = string(loss)
        flush(stderr)
        expected = ForwardDiff.gradient(loss, sample)
        @testset "Enzyme and Mooncake" begin
            for backend in (AutoEnzyme(; function_annotation = Enzyme.Const), AutoMooncake(; config = nothing))
                @test DifferentiationInterface.gradient(loss, backend, sample) ≈ expected rtol = 1.0e-7 atol = 1.0e-8
            end
        end
        @testset "Reactant CPU values and gradients" begin
            objective(x) = Reactant.@allowscalar loss(x)
            tracked = Reactant.to_rarray(sample)
            executable = Reactant.@compile objective(tracked)
            gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
            @test Reactant.to_number(executable(tracked)) ≈ loss(sample) rtol = 1.0e-8
            actual = Array(only(gradient(Enzyme.Reverse, objective, tracked)))
            @test actual ≈ expected rtol = 1.0e-7 atol = 1.0e-8
        end
    end
end

@testset "GWNoise effective PSD composition" begin
    curve = TabulatedPSD([0.0, 128.0], [2.0, 2.0])
    noise = GaussianNoise(curve; sampling_frequency = 256.0, nsamples = 256, minimum_frequency = 20.0)
    metric = FrequencyDomainMetric(frequencies(noise), (:H1,), psd(noise); minimum_frequency = 20.0, nsamples = 256)
    @test isfinite(optimal_snr(noise_fd(noise), metric))
end

@testset "Series host boundary with AD array kernels" begin
    epoch0 = 1126259460.0
    signal = FrequencySeries(pipeline_signal([35.8, 1.2]); frequencies = PIPE_F, epoch = epoch0, channels = (:H1, :L1))
    data = FrequencySeries(PIPE_DATA; frequencies = PIPE_F, epoch = epoch0, channels = (:H1, :L1))
    check_same_grid(signal, data)
    a = GWMatchedFiltering.check_array(signal, PIPE_METRIC)
    b = GWMatchedFiltering.check_array(data, PIPE_METRIC)
    @test a === parent(signal)
    @test log_likelihood_ratio(signal, data, PIPE_METRIC) == log_likelihood_ratio(a, b, PIPE_METRIC)
    # Metadata is validated once above. Each backend sees only arrays and a fixed metric.
    loss(values) = log_likelihood_ratio(values .* a, b, PIPE_METRIC)
    sample = ones(size(a))
    expected = ForwardDiff.gradient(loss, sample)
    for backend in (AutoEnzyme(; function_annotation = Enzyme.Const), AutoMooncake(; config = nothing))
        @test DifferentiationInterface.gradient(loss, backend, sample) ≈ expected rtol = 1.0e-7
    end
    objective(values) = Reactant.@allowscalar loss(values)
    tracked = Reactant.to_rarray(sample)
    executable = Reactant.@compile objective(tracked)
    gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    @test Reactant.to_number(executable(tracked)) ≈ loss(sample) rtol = 1.0e-8
    @test Array(only(gradient(Enzyme.Reverse, objective, tracked))) ≈ expected rtol = 1.0e-7
    # Reconstruct an array result on the host, retaining the validated grid.
    scale_signal(values) = values .* a
    array_executable = Reactant.@compile scale_signal(tracked)
    restored = with_data(signal, Array(array_executable(tracked)))
    @test check_same_grid(restored, signal) === nothing
    @test parent(restored) ≈ a
    curve = TabulatedPSD([0.0, 128.0], [2.0, 2.0])
    noise = GaussianNoise(curve; sampling_frequency = 256.0, nsamples = 256, minimum_frequency = 20.0)
    draw = noise_fd(FrequencySeries, noise; epoch = epoch0, channels = (:H1,))
    metric = FrequencyDomainMetric(frequencies(noise), (:H1,), psd(noise); minimum_frequency = 20.0, nsamples = 256)
    @test optimal_snr(draw, metric) == optimal_snr(parent(draw), metric)
end
