using TestItems

@testsnippet Setup_acceptance_ad begin
    using GWSeries
    using Test
    using FrequencyDomainWaveforms
    using Mooncake: Mooncake

    include("ad_helpers.jl")
end

@testitem "Waveform parameter gradients (acceptance)" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_ad] begin
    test_waveform_gradients(
        Pair{String, Any}[
            "ForwardDiff" => AutoForwardDiff(),
            "Mooncake" => AutoMooncake(; config = nothing),
        ]
    )
end

@testitem "Mooncake boundaries and cache reuse" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_ad] begin
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        generator = G(; minimum_frequency = 20.0)
        function loss(x)
            p = (;
                mass_1 = x[1], mass_2 = x[2], a_1 = x[3], a_2 = x[4],
                luminosity_distance = 500.0, theta_jn = 0.4, phase = 0.2,
            )
            h = hphc_fd(AD_GRID, p, generator)
            return 1.0e22 * (sum(real, h[:, 1]) + sum(imag, h[:, 2]))
        end
        # Cross component ordering, spline segments, waveform regions and cutoff bins.
        cache = Mooncake.prepare_gradient_cache(loss, [36.0, 29.0, 0.3, 0.2])
        for x in (
                [36.0, 29.0, 0.3, 0.2], [30.0, 30.0, 0.6, 0.3],
                [29.0, 36.0, 0.2, 0.3], [180.0, 120.0, 0.3, 0.2],
            )
            expected = ForwardDiff.gradient(loss, x)
            _, tangents = Mooncake.value_and_gradient!!(cache, loss, x)
            @test all(isfinite, tangents[2])
            @test tangents[2] ≈ expected rtol = 1.0e-8
        end
    end
    isco_loss(x) = FrequencyDomainWaveforms.xas_isco(x[1])
    cache = Mooncake.prepare_gradient_cache(isco_loss, [0.0])
    for a in (-0.7, -1.0e-7, 0.0, 1.0e-7, 0.7)
        _, tangents = Mooncake.value_and_gradient!!(cache, isco_loss, [a])
        @test tangents[2][1] ≈ ForwardDiff.derivative(FrequencyDomainWaveforms.xas_isco, a) rtol = 1.0e-10
    end
end
