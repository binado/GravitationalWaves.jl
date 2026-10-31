using TestItems

@testsnippet Setup_acceptance_ad begin
    using Test
    using GWSeries
    using GWDetectors
    using ForwardDiff
    using DifferentiationInterface
    using Mooncake: Mooncake

    include("ad_helpers.jl")
end

@testitem "GWDetectors AD integration" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_ad] begin
    samples = (AD_INPUT, AD_INPUT .+ [0.1, 0.02, -0.03, 0.125, 0.25, 0.01, -0.02, 0.04, -0.01, 0.001, 0.002, 10.0])
    @testset "Mooncake" begin
        for loss in (geometry_loss, antenna_loss, timing_loss, projection_loss)
            for sample in samples
                g = DifferentiationInterface.gradient(loss, AutoMooncake(; config = nothing), sample)
                @test g ≈ ForwardDiff.gradient(loss, sample) rtol = 1.0e-9 atol = 1.0e-10
            end
        end
    end
end
