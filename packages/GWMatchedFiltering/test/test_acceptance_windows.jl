using TestItems

@testsnippet Setup_acceptance_windows begin
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

    include("prototypes.jl")
end

@testitem "Fixed Planck window ForwardDiff and Mooncake" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_windows] begin
    td = time_prototype(; epoch = 1126259460.125, nsamples = 17, sampling_frequency = 8, channels = (:H1,))
    window = PlanckWindow(td; start_width = 0.5, end_width = 0.75)
    w = window.weights
    loss(h) = sum(abs2, apply_window(window, h))
    inplace_loss(h) = sum(abs2, apply_window!(copy(h), window))
    buffered_loss(h) = sum(abs2, apply_window!(similar(h), window, h))
    for h in (collect(range(-2.0, 3.0; length = 17)), reshape(collect(1.0:34), 17, 2), reshape(collect(1.0:68), 17, 2, 2))
        expected_value = sum(w[i[1]]^2 * h[i]^2 for i in CartesianIndices(h))
        expected_gradient = 2 .* w .^ 2 .* h
        for objective in (loss, inplace_loss, buffered_loss), backend in (AutoForwardDiff(), AutoMooncake(; config = nothing))
            value, gradient = DifferentiationInterface.value_and_gradient(objective, backend, h)
            @test value ≈ expected_value
            @test gradient ≈ expected_gradient
            @test size(gradient) == size(h)
            @test all(iszero, selectdim(gradient, 1, 1)) && all(iszero, selectdim(gradient, 1, length(w)))
        end
    end
end
