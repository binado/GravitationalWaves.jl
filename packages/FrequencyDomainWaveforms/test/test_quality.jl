using TestItems

@testsnippet Setup_quality begin
    using AbstractWaveformsTests

    using Aqua
    using JET
    import AbstractWaveforms
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(FrequencyDomainWaveforms)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    # JET's analysis is tied to the Julia version it was built for.
    if VERSION.major == 1 && VERSION.minor in (12, 13) && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(FrequencyDomainWaveforms; target_modules = (FrequencyDomainWaveforms,))
        end

        # The waveform kernels must be free of runtime dispatch for the AD backends.
        @testset "JET optimisation" begin
            f = range(0, 128; step = 0.5)
            params = (
                mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
                a_1 = 0.5, a_2 = 0.3, tilt_1 = 0.0, tilt_2 = Float64(π),
            )
            for generator in (TaylorF2(), IMRPhenomD(), IMRPhenomXAS())
                JET.@test_opt target_modules = (FrequencyDomainWaveforms, AbstractWaveforms) hphc_fd(f, params, generator)
                JET.@test_opt target_modules = (FrequencyDomainWaveforms,) amplitude_phase_fd(f, params, generator)
            end
        end
    end
end

@testitem "type stability" setup = [Common, Setup_quality] begin
    f = range(0, 128; step = 0.5)
    params = (
        chirp_mass = 28.0, mass_ratio = 0.8, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
    )
    for generator in (TaylorF2(), IMRPhenomD(), IMRPhenomXAS())
        @test (@inferred hphc_fd(f, params, generator)) isa Matrix
    end
end
