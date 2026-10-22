if !isdefined(@__MODULE__, :time_prototype)
    include(joinpath(@__DIR__, "../../packages/GWSeries/test/prototypes.jl"))
end
using Test
using GWMatchedFiltering
using GWSeries
using GWDetectors
using FrequencyDomainWaveforms
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
const PIPE_H1_DOMAIN = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:H1,))
const PIPE_L1_DOMAIN = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:L1,))
const PIPE_GENERATOR = TaylorF2(; minimum_frequency = 20.0, maximum_frequency = 128.0)
const PIPE_SKY = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
const PIPE_METRIC_domain = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:H1, :L1))
const PIPE_METRIC = FrequencyDomainMetric(PIPE_METRIC_domain, repeat(reshape([2.0 0.2; 0.2 1.0] .* 1.0e-46, 1, 2, 2), length(PIPE_F), 1, 1); minimum_frequency = 20.0)

function pipeline_signal(x)
    params = (; mass_1 = x[1], mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = x[2])
    pol = hphc_fd(PIPE_H1_DOMAIN, params, PIPE_GENERATOR)
    h1 = project_fd(pol, PIPE_SKY, detectors.H1; epoch = epoch(PIPE_H1_DOMAIN), channel = :H1)
    l1 = project_fd(pol, PIPE_SKY, detectors.L1; epoch = epoch(PIPE_L1_DOMAIN), channel = :L1)
    return hcat(parent(h1), parent(l1))
end

const PIPE_DATA = pipeline_signal([36.0, 1.3])
pipeline_loss(x) = log_likelihood_ratio(pipeline_signal(x), PIPE_DATA, PIPE_METRIC)

@testset "GWMatchedFiltering AD integration" begin
    for (loss, sample) in ((diagonal_loss, [0.8, 0.2]), (correlated_loss, [0.8, 0.2]), (statistics_loss, [0.8, 0.2]))
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

@testset "Waveform/projection/likelihood ForwardDiff and Mooncake" begin
    sample = [35.8, 1.2]
    expected = ForwardDiff.gradient(pipeline_loss, sample)
    actual = DifferentiationInterface.gradient(pipeline_loss, AutoMooncake(; config = nothing), sample)
    @test actual ≈ expected rtol = 1.0e-7 atol = 1.0e-8
end

@testset "GWNoise effective PSD composition" begin
    curve = TabulatedPSD([0.0, 128.0], [2.0, 2.0])
    noise = GaussianNoise(curve, time_prototype(; epoch = 0.0, nsamples = 256, sampling_frequency = 256.0, channels = (:H1,)); minimum_frequency = 20.0)
    metric_domain = frequency_prototype(; frequencies = frequencies(noise), epoch = 0.0, channels = (:H1,))
    metric = FrequencyDomainMetric(metric_domain, psd(noise); minimum_frequency = 20.0, nsamples = 256)
    @test isfinite(optimal_snr(noise_fd(noise), metric))
end


@testset "Series metadata at the inference boundary" begin
    fd = PIPE_METRIC_domain
    a = pipeline_signal([35.8, 1.2]);b = PIPE_DATA
    check_data(fd, a)
    check_data(fd, b)
    check_compatible(fd, PIPE_METRIC)
    check_template(a, PIPE_METRIC)
    @test check_same_grid(fd, frequency_prototype(; frequencies = PIPE_F, epoch = epoch(fd), channels = channels(fd))) === nothing
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
end

include("test_vectors.jl")
include("test_windows.jl")

@testset "Flatten, raw kernel, and series reconstruction AD" begin
    data = FrequencySeries(PIPE_DATA; metadata(PIPE_METRIC_domain)...)
    function boundary_loss(values)
        array, meta = flatten(data)
        signal = FrequencySeries(values .* array; meta...)
        return log_likelihood_ratio(parent(signal), array, PIPE_METRIC)
    end
    # Stay away from the likelihood maximum, where the exact gradient is zero
    # and a relative-only comparison measures cancellation roundoff.
    sample = fill(0.8, size(data))
    expected = ForwardDiff.gradient(boundary_loss, sample)
    for backend in (AutoEnzyme(; function_annotation = Enzyme.Const), AutoMooncake(; config = nothing))
        @test DifferentiationInterface.gradient(boundary_loss, backend, sample) ≈ expected rtol = 1.0e-7
    end
    objective(values) = Reactant.@allowscalar boundary_loss(values)
    tracked = Reactant.to_rarray(sample)
    executable = Reactant.@compile objective(tracked)
    gradient = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    @test Reactant.to_number(executable(tracked)) ≈ boundary_loss(sample) rtol = 1.0e-8
    @test Array(only(gradient(Enzyme.Reverse, objective, tracked))) ≈ expected rtol = 1.0e-7
end
