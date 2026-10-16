using ForwardDiff

@testset "spins" begin
    @testset "effective spins" begin
        @test effective_spin(0.5, 0.5, 0.7) ≈ 0.5
        @test effective_spin(0.4, -0.4, 1.0) == 0.0
        @test effective_spin(0.3, 0.2, 0.5) ≈ (0.3 + 0.2 * 0.5) / 1.5
        # Spin 1 dominates when q = 1 and the spins are equal.
        @test effective_precessing_spin(0.3, 0.4, 0.0, 0.0, 0.8) ≈ 0.5
        @test effective_precessing_spin(0.0, 0.0, 0.6, 0.8, 1.0) ≈ 1.0
        @test effective_precessing_spin(0.0, 0.0, 0.0, 0.0, 0.5) == 0.0
    end

    @testset "alignment test uses the rounded double" begin
        @test is_aligned(0.0, 1.0)
        @test is_aligned(0.5, 0.0)
        @test is_aligned(0.5, Float64(π))
        @test !is_aligned(0.5, 1.0)
        @test !is_aligned(0.5, π - 1.0e-12)
    end
end

@testset "frames" begin
    base = (mass_1 = 36.0, mass_2 = 29.0, theta_jn = 0.4, phase = 1.3)
    precessing = merge(
        base, (a_1 = 0.6, a_2 = 0.4, tilt_1 = 0.7, tilt_2 = 1.9, phi_12 = 1.1, phi_jl = 2.5)
    )

    @testset "aligned branch is exact" begin
        s = spins_to_lal_frame(merge(base, (a_1 = 0.3, a_2 = 0.2, tilt_2 = Float64(π))), 50.0)
        @test s.iota == 0.4
        @test (s.s1x, s.s1y, s.s2x, s.s2y) == (0.0, 0.0, 0.0, 0.0)
        @test s.s1z == 0.3
        @test s.s2z == 0.2 * cos(Float64(π)) == -0.2
    end

    @testset "zero spin magnitude counts as aligned" begin
        s = spins_to_lal_frame(merge(base, (a_1 = 0.0, a_2 = 0.5, tilt_1 = 1.0)), 50.0)
        @test s.s1z == 0.0
        @test s.iota == 0.4
    end

    @testset "precessing invariants" begin
        s = spins_to_lal_frame(precessing, 50.0)
        p = precessing
        @test hypot(s.s1x, s.s1y, s.s1z) ≈ p.a_1
        @test hypot(s.s2x, s.s2y, s.s2z) ≈ p.a_2
        # The tilts are angles to L̂, which is the z axis of the LAL frame.
        @test s.s1z ≈ p.a_1 * cos(p.tilt_1)
        @test s.s2z ≈ p.a_2 * cos(p.tilt_2)
        # The relative azimuth of the spins is phi_12.
        s1s2 = s.s1x * s.s2x + s.s1y * s.s2y + s.s1z * s.s2z
        @test s1s2 ≈ p.a_1 * p.a_2 * (sin(p.tilt_1) * sin(p.tilt_2) * cos(p.phi_12) + cos(p.tilt_1) * cos(p.tilt_2))
        @test 0 <= s.iota <= π
    end

    @testset "zero spins give iota = theta_jn" begin
        # Nearly aligned spins: J ≈ L so the inclination is theta_jn.
        s = transform_precessing_initial_conditions(0.4, 1.0, 1.0e-9, 1.0e-9, 0.5, 0.5, 0.5, 36.0, 29.0, 50.0, 1.3)
        @test s.iota ≈ 0.4 atol = 1.0e-6
    end

    @testset "phase rotates the spins about L" begin
        a = spins_to_lal_frame(precessing, 50.0)
        b = spins_to_lal_frame(merge(precessing, (phase = precessing.phase + 0.5,)), 50.0)
        @test a.iota ≈ b.iota
        @test a.s1z ≈ b.s1z
        @test hypot(a.s1x, a.s1y) ≈ hypot(b.s1x, b.s1y)
        # N̂ sits at azimuth π/2 - phase, so a phase shift of δ rotates the in-plane spins by -δ about z.
        @test mod2pi(atan(b.s1y, b.s1x) - atan(a.s1y, a.s1x) + π) - π ≈ -0.5
    end

    @testset "N parallel to L falls back to a fixed azimuth reference" begin
        s = transform_precessing_initial_conditions(0.0, 1.0, 0.7, 1.9, 1.1, 0.6, 0.4, 36.0, 29.0, 50.0, 1.3)
        @test all(isfinite, values(s))
        @test hypot(s.s1x, s.s1y, s.s1z) ≈ 0.6
    end

    @testset "automatic differentiation" begin
        fd(f, x; h = 1.0e-6) = (f(x + h) - f(x - h)) / 2h
        for key in (:iota, :s1x, :s1y, :s1z, :s2x, :s2y, :s2z)
            f(tilt) = getproperty(spins_to_lal_frame(merge(precessing, (tilt_1 = tilt,)), 50.0), key)
            @test ForwardDiff.derivative(f, 0.7) ≈ fd(f, 0.7) rtol = 1.0e-6 atol = 1.0e-8
            g(phi) = getproperty(spins_to_lal_frame(merge(precessing, (phi_jl = phi,)), 50.0), key)
            @test ForwardDiff.derivative(g, 2.5) ≈ fd(g, 2.5) rtol = 1.0e-6 atol = 1.0e-8
        end
        @test @inferred(spins_to_lal_frame(precessing, 50.0)) isa NamedTuple
    end
end
