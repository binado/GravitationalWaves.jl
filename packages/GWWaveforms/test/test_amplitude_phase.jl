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
            test_amplitude_phase_contract(g; params, f)
            a = @inferred amplitude_phase(f, params, g)
            @test eltype(a.amplitude) === eltype(a.phase) === Float64
            active = a.amplitude .> 0
            @test any(abs.(a.phase[active]) .> 2π)
            shifted = amplitude_phase(f, merge(params, (; phase = 0.4)), g)
            @test shifted.amplitude == a.amplitude
            @test all(isapprox.(shifted.phase[active] .- a.phase[active], 0.4; atol = 1.0e-10))
            distant = amplitude_phase(f, merge(params, (; luminosity_distance = 1000.0)), g)
            @test distant.amplitude == a.amplitude ./ 2
            @test distant.phase == a.phase
            empty = amplitude_phase(f, params, G(; minimum_frequency = 600.0, maximum_frequency = 700.0))
            @test all(iszero, empty.amplitude) && all(iszero, empty.phase)
            @test_throws ArgumentError amplitude_phase(1:1:100, params, g)
            @test_throws ArgumentError amplitude_phase(f, params, G(; maximum_frequency = 0.0))
        end
    end
end
