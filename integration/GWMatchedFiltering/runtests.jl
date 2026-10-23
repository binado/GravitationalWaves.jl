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
using Random
using ForwardDiff
using DifferentiationInterface
using Mooncake: Mooncake

include("../../packages/GWMatchedFiltering/test/whitening_helpers.jl")
include("../../packages/GWMatchedFiltering/test/ad_helpers.jl")

const PIPE_F = 0.0:1.0:128.0
const PIPE_H1_DOMAIN = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:H1,))
const PIPE_L1_DOMAIN = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:L1,))
const PIPE_GENERATOR = TaylorF2(; minimum_frequency = 20.0, maximum_frequency = 128.0)
const PIPE_SKY = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
const PIPE_METRIC_domain = frequency_prototype(; frequencies = PIPE_F, epoch = 1126259460.0, channels = (:H1, :L1))
const PIPE_DF = step(PIPE_F)
# Correlated H1/L1 noise above 20 Hz; zero blocks mark the inactive bins.
const PIPE_W = whitening_blocks(
    cat((f >= 20.0 ? ComplexF64[2.0 0.2; 0.2 1.0] .* 1.0e-46 : zeros(ComplexF64, 2, 2) for f in PIPE_F)...; dims = 3)
)

function pipeline_signal(x)
    params = (; mass_1 = x[1], mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = x[2])
    pol = hphc_fd(PIPE_H1_DOMAIN, params, PIPE_GENERATOR)
    h1 = project_fd(pol, PIPE_SKY, detectors.H1; epoch = epoch(PIPE_H1_DOMAIN), channel = :H1)
    l1 = project_fd(pol, PIPE_SKY, detectors.L1; epoch = epoch(PIPE_L1_DOMAIN), channel = :L1)
    return hcat(parent(h1), parent(l1))
end

const PIPE_DATA = pipeline_signal([36.0, 1.3])
pipeline_loss(x) = log_likelihood_ratio(pipeline_signal(x), PIPE_DATA, PIPE_W, PIPE_DF)

@testset "GWMatchedFiltering AD integration" begin
    for (loss, sample) in ((diagonal_loss, [0.8, 0.2]), (correlated_loss, [0.8, 0.2]), (statistics_loss, [0.8, 0.2]))
        @info "Checking GWMatchedFiltering derivatives" objective = string(loss)
        flush(stderr)
        expected = ForwardDiff.gradient(loss, sample)
        @testset "Mooncake" begin
            @test DifferentiationInterface.gradient(loss, AutoMooncake(; config = nothing), sample) ≈ expected rtol = 1.0e-7 atol = 1.0e-8
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
    f = range(0.0; step = 1.0, length = 129)
    S = reshape(band_psd(curve, f; fmin = 20.0), :, 1)
    W = whitening_factor(S)
    @test isfinite(optimal_snr(sample_noise(Xoshiro(1), S, step(f)), W, step(f)))
end


@testset "Series metadata at the inference boundary" begin
    fd = PIPE_METRIC_domain
    a = pipeline_signal([35.8, 1.2]);b = PIPE_DATA
    check_data(fd, a)
    check_data(fd, b)
    check_template(a, PIPE_W, PIPE_DF)
    @test check_same_grid(fd, frequency_prototype(; frequencies = PIPE_F, epoch = epoch(fd), channels = channels(fd))) === nothing
    loss(values) = log_likelihood_ratio(values .* a, b, PIPE_W, PIPE_DF)
    sample = ones(size(a))
    expected = ForwardDiff.gradient(loss, sample)
    @test DifferentiationInterface.gradient(loss, AutoMooncake(; config = nothing), sample) ≈ expected rtol = 1.0e-7
end

include("test_vectors.jl")
include("test_windows.jl")

@testset "Flatten, raw kernel, and series reconstruction AD" begin
    data = FrequencySeries(PIPE_DATA; metadata(PIPE_METRIC_domain)...)
    function boundary_loss(values)
        array, meta = flatten(data)
        signal = FrequencySeries(values .* array; meta...)
        return log_likelihood_ratio(parent(signal), array, PIPE_W, PIPE_DF)
    end
    # Stay away from the likelihood maximum, where the exact gradient is zero
    # and a relative-only comparison measures cancellation roundoff.
    sample = fill(0.8, size(data))
    expected = ForwardDiff.gradient(boundary_loss, sample)
    @test DifferentiationInterface.gradient(boundary_loss, AutoMooncake(; config = nothing), sample) ≈ expected rtol = 1.0e-7
end
