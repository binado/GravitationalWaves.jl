# Both backends are loaded unqualified: the interface names are shared, not ambiguous.
using LALWaveforms

@testset "common interface" begin
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
        a_1 = 0.3, a_2 = 0.1, tilt_1 = 0.0, tilt_2 = 0.0,
    )
    f = range(0, 512; step = 0.25)
    # Code written against the interface only.
    strain_norm(g) = sqrt(sum(abs2, hphc_fd(f, params, g)[:, 1]))

    if lal_available()
        @test AbstractWaveformGenerator === LALWaveforms.AbstractWaveformGenerator
        @test FrequencyDomainWaveforms.hphc_fd === LALWaveforms.hphc_fd
        @test FrequencyDomainWaveforms.validate_parameters === LALWaveforms.validate_parameters
        native, lal = IMRPhenomD(), LALWaveformGenerator("IMRPhenomD")
        @test approximant(native) == approximant(lal)
        @test native isa AbstractWaveformGenerator && lal isa AbstractWaveformGenerator
        @test strain_norm(native) ≈ strain_norm(lal) rtol = 1.0e-10
        @test is_differentiable(native) && !is_differentiable(lal)
    end
end
