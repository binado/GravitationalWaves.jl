@testset "generators and input validation" begin
    base = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
    )
    f = range(0, 512; step = 0.25)

    @testset "generators" begin
        g = IMRPhenomD()
        @test g isa AbstractWaveformGenerator
        @test minimum_frequency(g) == 20.0
        @test reference_frequency(g) == 50.0
        @test maximum_frequency(g) === nothing
        @test approximant(g) == "IMRPhenomD"
        @test approximant(TaylorF2()) == "TaylorF2"
        @test spin_support(g) === :aligned
        @test is_differentiable(g) && !supports_mode_array(g)
        @test TaylorF2(; reference_frequency = 20).reference_frequency === 20.0
        @test_throws ArgumentError TaylorF2(; pn_amplitude_order = 1)
        @test_throws ArgumentError TaylorF2(; pn_phase_order = 8)
        @test_throws ArgumentError TaylorF2(; pn_spin_order = -2)
        @test_throws ArgumentError IMRPhenomD(; pn_amplitude_order = 1)
        # Options a model does not use are not fields of its generator.
        @test_throws MethodError IMRPhenomD(; mode_array = [(2, 2)])
        @test_throws MethodError IMRPhenomD(; pn_phase_order = 3)
        @test_throws MethodError TaylorF2(; pn_tidal_order = 3)
        @test sprint(show, TaylorF2(; maximum_frequency = 300)) ==
            "TaylorF2(\"TaylorF2\"; reference_frequency=50.0, minimum_frequency=20.0, maximum_frequency=300.0)"
    end

    @testset "parameters" begin
        s = IMRPhenomD()
        @test validate_parameters(f, base, s) === nothing
        @test_throws ArgumentError validate_parameters(f, merge(base, (; foo = 1.0)), s)
        @test_throws ArgumentError validate_parameters(f, Base.structdiff(base, (; phase = 0.0)), s)
        # Precessing spins and eccentricity are outside both models.
        @test_throws ArgumentError validate_parameters(f, merge(base, (; a_1 = 0.5, tilt_1 = 0.7)), s)
        @test_throws ArgumentError validate_parameters(f, merge(base, (; eccentricity = 0.1)), s)
        @test validate_parameters(f, merge(base, (; a_1 = 0.5, tilt_1 = Float64(π))), s) === nothing
        # A spin of zero is aligned whatever its tilt.
        @test hphc_fd(f, merge(base, (; a_1 = 0.0, tilt_1 = 0.7)), s).plus ==
            hphc_fd(f, base, s).plus
        # Other mass parameterisations convert.
        q = base.mass_2 / base.mass_1
        mc = (base.mass_1 * base.mass_2)^(3 / 5) / (base.mass_1 + base.mass_2)^(1 / 5)
        alt = (;
            chirp_mass = mc, mass_ratio = q, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
        )
        @test hphc_fd(f, alt, s).plus ≈ hphc_fd(f, base, s).plus rtol = 1.0e-12
    end

    @testset "parameter values" begin
        for G in (IMRPhenomD, TaylorF2, IMRPhenomXAS)
            g = G(; minimum_frequency = 30.0, maximum_frequency = 400.0)
            @test validate_parameters(f, base, g) === nothing
            @test_throws ArgumentError validate_parameters(f, merge(base, (; mass_2 = -1.0)), g)
            @test_throws ArgumentError validate_parameters(f, merge(base, (; luminosity_distance = 0.0)), g)
            if g isa TaylorF2
                @test_throws ArgumentError validate_parameters(range(1, 512; step = 0.25), base, g)
            else
                @test validate_parameters(range(1, 512; step = 0.25), base, g) === nothing
            end
            # Parameter checks are explicit through `validate_parameters`.
            @test_throws ArgumentError validate_parameters(f, base, G(; minimum_frequency = 0.0))
            @test_throws ArgumentError validate_parameters(f, base, G(; minimum_frequency = 450.0, maximum_frequency = 400.0))
            @test_throws ArgumentError validate_parameters(f, base, G(; reference_frequency = -1.0))
        end
        g = IMRPhenomD()
        for G in (IMRPhenomD, IMRPhenomXAS)
            @test_throws ArgumentError validate_parameters(f, merge(base, (; a_1 = 1.5)), G())
        end
        @test_throws ArgumentError validate_parameters(f, merge(base, (; mass_1 = 5000.0, mass_2 = 5000.0)), g)
        # The default cutoff depends on the masses, which the static path does not support.
        for z in (TaylorF2(; maximum_frequency = 0), IMRPhenomD(; maximum_frequency = 0))
            @test_throws ArgumentError validate_parameters(f, base, z)
        end
    end

    @testset "grid" begin
        s = IMRPhenomD()
        @test validate_parameters(range(1, 512; step = 0.25), base, s) === nothing
        @test_throws ArgumentError validate_parameters(range(0, 512; step = -0.25), base, s)
        @test_throws ArgumentError validate_parameters(range(0, 0; length = 1), base, s)
    end

    @testset "masking and f_max" begin
        for G in (IMRPhenomD, TaylorF2, IMRPhenomXAS)
            h = hphc_fd(f, base, G(; minimum_frequency = 30, maximum_frequency = 200))
            below = f .< 30
            above = f .> 200
            inside = .!(below .| above)
            @test all(iszero, h.plus[below]) && all(iszero, h.cross[below])
            @test all(iszero, h.plus[above]) && all(iszero, h.cross[above])
            @test all(!iszero, h.plus[inside][1:(end - 1)])
            # `maximum_frequency = nothing` means the last frequency of the grid.
            h_default = hphc_fd(f, base, G(; minimum_frequency = 30))
            h_explicit = hphc_fd(
                f, base, G(; minimum_frequency = 30, maximum_frequency = last(f))
            )
            @test h_default.plus == h_explicit.plus
            @test h_default.cross == h_explicit.cross
        end
    end

    @testset "orientation" begin
        # Face-on: hc = -i hp; edge-on: no cross polarisation, hp halved.
        for G in (IMRPhenomD, TaylorF2, IMRPhenomXAS)
            s = G()
            face = hphc_fd(f, merge(base, (; theta_jn = 0.0)), s)
            @test face.cross ≈ -im .* face.plus
            edge = hphc_fd(f, merge(base, (; theta_jn = π / 2)), s)
            @test maximum(abs, edge.cross) < 1.0e-12 * maximum(abs, edge.plus)
            @test edge.plus ≈ face.plus ./ 2
        end
    end

    @testset "element type" begin
        s = TaylorF2()
        @test eltype(hphc_fd(f, base, s).plus) === ComplexF64
        @test eltype(hphc_fd(f, merge(base, (; mass_1 = 36)), s).plus) === ComplexF64
    end
end
