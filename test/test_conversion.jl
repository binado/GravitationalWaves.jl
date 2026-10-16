@testset "conversion" begin
    LW = LALWaveforms
    base = (mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3)

    @testset "aligned branch is exact" begin
        p = LW.BBHParameters(merge(base, (a_1 = 0.3, a_2 = 0.2, tilt_2 = Float64(π))))
        s = LW.bilby_to_lalsimulation_spins(p, 50.0)
        @test s.iota == 0.4
        @test (s.s1x, s.s1y, s.s2x, s.s2y) == (0.0, 0.0, 0.0, 0.0)
        @test s.s1z == 0.3
        @test s.s2z == 0.2 * cos(Float64(π)) == -0.2
    end

    @testset "zero spin magnitude counts as aligned" begin
        p = LW.BBHParameters(merge(base, (a_1 = 0.0, a_2 = 0.5, tilt_1 = 1.0)))
        @test LW.bilby_to_lalsimulation_spins(p, 50.0).s1z == 0.0
    end

    @testset "LAL branch agrees with aligned for non-precessing input" begin
        p = LW.BBHParameters(merge(base, (a_1 = 0.3, a_2 = 0.2, tilt_2 = Float64(π))))
        a = LW.bilby_to_lalsimulation_spins(p, 50.0)
        b = LW.transform_precessing(p, 50.0)
        for k in keys(a)
            @test getproperty(a, k) ≈ getproperty(b, k) atol = 1.0e-12
        end
    end

    @testset "precessing input goes through LAL" begin
        p = LW.BBHParameters(
            merge(base, (a_1 = 0.6, a_2 = 0.4, tilt_1 = 0.7, tilt_2 = 1.9, phi_12 = 1.1, phi_jl = 2.5))
        )
        s = LW.bilby_to_lalsimulation_spins(p, 50.0)
        @test hypot(s.s1x, s.s1y, s.s1z) ≈ 0.6
        @test hypot(s.s2x, s.s2y, s.s2z) ≈ 0.4
    end
end
