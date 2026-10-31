using TestItems

@testsnippet Setup_interface begin
    using AbstractWaveformsTests

    # Both backends are loaded unqualified: the interface names are shared, not ambiguous.
    using LALWaveforms

    # The numeric cross-check needs the LAL shared libraries (`LALWAVEFORMS_LIBDIR`); the rest of the
    # interface checks do not. The full LAL comparison is `just test-reference-waveforms`.
    lal_available() = try
        hphc_fd(
            range(0, 64; step = 1.0),
            (mass_1 = 30.0, mass_2 = 30.0, luminosity_distance = 100.0, theta_jn = 0.0, phase = 0.0),
            LALWaveformGenerator("TaylorF2"; minimum_frequency = 20.0),
        )
        true
    catch err
        err isa InterruptException && rethrow()
        false
    end
end

@testitem "common interface" setup = [Common, Setup_interface] begin
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
        a_1 = 0.3, a_2 = 0.1, tilt_1 = 0.0, tilt_2 = 0.0,
    )
    f = range(0, 512; step = 0.25)
    # Code written against the interface only.
    strain_norm(g) = sqrt(sum(abs2, hphc_fd(f, params, g)[:, 1]))

    @test AbstractWaveformGenerator === LALWaveforms.AbstractWaveformGenerator
    @test FrequencyDomainWaveforms.hphc_fd === LALWaveforms.hphc_fd
    @test FrequencyDomainWaveforms.validate_parameters === LALWaveforms.validate_parameters
    if lal_available()
        native, lal = IMRPhenomD(), LALWaveformGenerator("IMRPhenomD")
        @test approximant(native) == approximant(lal)
        @test native isa AbstractWaveformGenerator && lal isa AbstractWaveformGenerator
        @test strain_norm(native) ≈ strain_norm(lal) rtol = 1.0e-10
        @test is_differentiable(native) && !is_differentiable(lal)
    end
end
