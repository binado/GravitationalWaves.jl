using ForwardDiff

to_masses(nt) = convert_parameters(nt, (:mass_1, :mass_2))

@testset "convert_parameters" begin
    m1, m2 = 36.0, 29.0
    mc, mt, q, eta = chirp_mass(m1, m2), total_mass(m1, m2), mass_ratio(m1, m2), symmetric_mass_ratio(m1, m2)
    masses = (:mass_1, :mass_2)

    @testset "routes to component masses" begin
        for nt in (
                (chirp_mass = mc, mass_ratio = q), (total_mass = mt, mass_ratio = q),
                (chirp_mass = mc, total_mass = mt), (symmetric_mass_ratio = eta, total_mass = mt),
                (mass_1 = m1, chirp_mass = mc), (mass_2 = m2, chirp_mass = mc),
            )
            out = convert_parameters(nt, masses)
            @test out.mass_1 ≈ m1
            @test out.mass_2 ≈ m2
        end
    end

    @testset "only the needed conversions run" begin
        out = convert_parameters((chirp_mass = mc, mass_ratio = q), masses)
        @test keys(out) == (:chirp_mass, :mass_ratio, :mass_1, :mass_2)
        out = convert_parameters((mass_1 = m1, mass_2 = m2), (:chirp_mass,))
        @test keys(out) == (:mass_1, :mass_2, :chirp_mass)
    end

    @testset "given keys are kept and extras pass through" begin
        nt = (mass_1 = m1, mass_2 = m2, theta_jn = 0.4)
        @test convert_parameters(nt, masses) === nt
        out = convert_parameters(nt, (:mass_ratio,))
        @test out.theta_jn == 0.4
        @test out.mass_ratio == m2 / m1
    end

    @testset "source frame" begin
        out = convert_parameters((mass_1_source = 30.0, mass_2_source = 20.0, redshift = 0.5), masses)
        @test out.mass_1 == 45.0
        @test out.mass_2 == 30.0
        out = convert_parameters((mass_1 = 45.0, mass_2 = 30.0, redshift = 0.5), (:mass_1_source, :mass_2_source))
        @test out.mass_1_source ≈ 30.0
        # Chained: chirp mass and mass ratio in the source frame, to detector component masses.
        src = (chirp_mass_source = mc / 1.5, mass_ratio = q, redshift = 0.5)
        out = convert_parameters(src, masses)
        @test out.mass_1 ≈ m1
    end

    @testset "spins" begin
        nt = (mass_1 = m1, mass_2 = m2, a_1 = 0.6, a_2 = 0.4, tilt_1 = 0.7, tilt_2 = 1.9)
        out = convert_parameters(nt, (:chi_eff, :chi_p))
        @test out.chi_eff ≈ effective_spin(0.6 * cos(0.7), 0.4 * cos(1.9), q)
        @test out.chi_p ≈ effective_precessing_spin(0.6 * sin(0.7), 0.4 * sin(1.9), q)
    end

    @testset "redundant inputs are checked" begin
        @test convert_parameters((mass_1 = m1, mass_2 = m2, chirp_mass = mc, mass_ratio = q), masses).mass_1 == m1
        err = try
            convert_parameters((mass_1 = m1, mass_2 = m2, chirp_mass = 1.1 * mc), masses)
        catch e
            e
        end
        @test err isa InconsistentParametersError
        @test err.key == :chirp_mass
        @test err.from == (:mass_1, :mass_2)
        @test occursin("chirp_mass", sprint(showerror, err))
        # A loose enough tolerance accepts the same input.
        @test convert_parameters((mass_1 = m1, mass_2 = m2, chirp_mass = 1.1 * mc), masses; rtol = 0.5).mass_1 == m1
        # A zero that the other keys imply to within rounding is not flagged.
        @test convert_parameters((a_1 = 0.5, tilt_1 = π / 2, spin_1z = 0.0), (:spin_1z,)).spin_1z == 0.0
    end

    @testset "underdetermined targets" begin
        err = try
            convert_parameters((chirp_mass = mc,), masses)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("mass_1", err.msg) && occursin("chirp_mass", err.msg)
        @test_throws ArgumentError convert_parameters((mass_1 = m1,), (:redshift,))
    end

    @testset "type stability and AD" begin
        # The plan resolves at compile time only if `targets` is a literal in the compiled code, so
        # check through a wrapper (@inferred evaluates its arguments first).
        @test @inferred(to_masses((chirp_mass = mc, mass_ratio = q))) isa NamedTuple
        f(x) = convert_parameters((chirp_mass = x, mass_ratio = q), masses).mass_1
        h = 1.0e-6
        @test ForwardDiff.derivative(f, mc) ≈ (f(mc + h) - f(mc - h)) / 2h rtol = 1.0e-6
        d(x) = convert_parameters((mass_1 = x, mass_2 = m2), (:chirp_mass,)).chirp_mass
        @test ForwardDiff.derivative(d, m1) ≈ ForwardDiff.derivative(m -> chirp_mass(m, m2), m1)
    end
end
