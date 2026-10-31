using TestItems

@testsnippet Setup_masses begin
    using ForwardDiff
end

@testitem "masses" setup = [Common, Setup_masses] begin
    G = GWParameterConversion
    pairs = [(36.0, 29.0), (50.0, 20.0), (30.0, 30.0), (10.5, 1.4), (1.4, 1.35)]

    @testset "closed forms" begin
        @test chirp_mass(30.0, 30.0) ≈ 30 * 2^(-1 / 5)
        @test total_mass(36.0, 29.0) == 65.0
        @test mass_ratio(40.0, 10.0) == 0.25
        @test symmetric_mass_ratio(30.0, 30.0) == 0.25
        @test reduced_mass(40.0, 10.0) == 8.0
        @test symmetric_mass_ratio(1.0, 1.0 + 1.0e-17) <= 0.25
    end

    @testset "round trips" begin
        for (m1, m2) in pairs
            mc, mt, q, eta = chirp_mass(m1, m2), total_mass(m1, m2), mass_ratio(m1, m2), symmetric_mass_ratio(m1, m2)
            @test all(
                total_mass_and_mass_ratio_to_component_masses(mt, q) .≈ (m1, m2)
            )
            @test all(
                chirp_mass_and_mass_ratio_to_component_masses(mc, q) .≈ (m1, m2)
            )
            @test chirp_mass_and_mass_ratio_to_total_mass(mc, q) ≈ mt
            @test chirp_mass_and_total_mass_to_symmetric_mass_ratio(mc, mt) ≈ eta
            @test symmetric_mass_ratio_to_mass_ratio(eta) ≈ q
            @test mass_1_and_chirp_mass_to_mass_ratio(m1, mc) ≈ q
            @test 0 < q <= 1 && 0 < eta <= 0.25
        end
    end

    @testset "mass_2 solve" begin
        @test mass_2_and_chirp_mass_to_mass_ratio(29.0, chirp_mass(36.0, 29.0)) ≈ mass_ratio(36.0, 29.0)
        # No real root in bilby's closed form once (chirp_mass / mass_2)^5 > 27/4.
        @test_throws DomainError mass_2_and_chirp_mass_to_mass_ratio(1.4, chirp_mass(10.5, 1.4))
    end

    @testset "redshift" begin
        @test source_to_detector_mass(30.0, 0.5) == 45.0
        @test detector_to_source_mass(source_to_detector_mass(30.0, 0.3), 0.3) ≈ 30.0
    end

    @testset "automatic differentiation" begin
        d(f, x) = ForwardDiff.derivative(f, x)
        # ∂Mc/∂m1 at fixed m2, against a central difference.
        h = 1.0e-6
        fd = (chirp_mass(36.0 + h, 29.0) - chirp_mass(36.0 - h, 29.0)) / 2h
        @test d(m -> chirp_mass(m, 29.0), 36.0) ≈ fd rtol = 1.0e-8
        for f in (
                m -> total_mass(m, 29.0), m -> mass_ratio(m, 29.0), m -> symmetric_mass_ratio(m, 29.0),
                m -> reduced_mass(m, 29.0), m -> total_mass_and_mass_ratio_to_component_masses(m, 0.8)[2],
                m -> chirp_mass_and_mass_ratio_to_component_masses(m, 0.8)[1],
                m -> chirp_mass_and_mass_ratio_to_total_mass(m, 0.8),
                m -> chirp_mass_and_total_mass_to_symmetric_mass_ratio(20.0, m + 30),
                m -> symmetric_mass_ratio_to_mass_ratio(m / 100),
                m -> mass_1_and_chirp_mass_to_mass_ratio(36.0, m),
                m -> mass_2_and_chirp_mass_to_mass_ratio(m + 29, 28.0),
                m -> source_to_detector_mass(m, 0.5), m -> detector_to_source_mass(m, 0.5),
            )
            @test isfinite(d(f, 20.0))
        end
        @test @inferred(chirp_mass(ForwardDiff.Dual(36.0, 1.0), 29.0)) isa ForwardDiff.Dual
    end
end
