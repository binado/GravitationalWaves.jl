@testset "analytic amplitude and phase" begin
    f = range(0, 512; step = 0.5)
    params = (;
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
        theta_jn = 0.4, phase = 0.2, a_1 = 0.5, a_2 = 0.3,
    )
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        @testset "$G" begin
            @test G <: AlignedSpinQuadrupolarWaveformGenerator
            g = G(; minimum_frequency = 20.1, maximum_frequency = 400.3)
            test_amplitude_phase_fd_contract(g; params, f)
            a = @inferred amplitude_phase_fd(f, params, g)
            @test eltype(a.amplitude) === eltype(a.phase) === Float64
            active = a.amplitude .> 0
            @test any(abs.(a.phase[active]) .> 2π)
            shifted = amplitude_phase_fd(f, merge(params, (; phase = 0.4)), g)
            @test shifted.amplitude == a.amplitude
            @test all(isapprox.(shifted.phase[active] .- a.phase[active], 0.4; atol = 1.0e-10))
            distant = amplitude_phase_fd(f, merge(params, (; luminosity_distance = 1000.0)), g)
            @test distant.amplitude == a.amplitude ./ 2
            @test distant.phase == a.phase
            empty = amplitude_phase_fd(f, params, G(; minimum_frequency = 600.0, maximum_frequency = 700.0))
            @test all(iszero, empty.amplitude) && all(iszero, empty.phase)
            @test_throws ArgumentError validate_parameters(f, params, G(; maximum_frequency = 0.0))
        end
    end
end

@testset "Phenom frequency ranges starting above zero" begin
    params = (;
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0,
        theta_jn = 0.4, phase = 0.2, a_1 = 0.5, a_2 = 0.3,
    )
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        g = G(; minimum_frequency = 20.1, maximum_frequency = 400.3)
        full = 0.0:0.125:512.0
        cropped = full[164:end]
        @test validate_parameters(cropped, params, g) === nothing
        a = amplitude_phase_fd(full, params, g)
        b = @inferred amplitude_phase_fd(cropped, params, g)
        @test b.amplitude == a.amplitude[164:end]
        @test b.phase == a.phase[164:end]
        h = hphc_fd(full, params, g)
        cropped_h = hphc_fd(cropped, params, g)
        @test cropped_h.plus == h.plus[164:end]
        @test cropped_h.cross == h.cross[164:end]

        # Offset by a fraction of the spacing: these frequencies are not integer bins of Δf.
        # TaylorF2's phase carries LAL's `-1/Δf` time shift, so it needs the same step as `fine`.
        offset = G === TaylorF2 ? (20.0625:0.0625:399.9375) : (20.0625:0.125:399.9375)
        fine = 0.0:0.0625:512.0
        expected = amplitude_phase_fd(fine, params, g)
        actual = amplitude_phase_fd(offset, params, g)
        indices = [searchsortedfirst(fine, fi) for fi in offset]
        @test actual.amplitude ≈ expected.amplitude[indices] rtol = 1.0e-12
        @test actual.phase ≈ expected.phase[indices] rtol = 1.0e-12
        for bad in (range(-1, 512; step = 0.25), 100.0:-1.0:20.0, range(20.0; step = 1.0, length = 1))
            @test_throws ArgumentError validate_parameters(bad, params, g)
        end
    end
end
