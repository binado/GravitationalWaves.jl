@testset "settings and input validation" begin
    base = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
    )
    f = range(0, 512; step = 0.25)

    @testset "WaveformSettings" begin
        s = WaveformSettings("IMRPhenomD")
        @test s.minimum_frequency == 20.0
        @test s.reference_frequency == 50.0
        @test s.maximum_frequency === nothing
        @test WaveformSettings(:TaylorF2) isa WaveformSettings{GWWaveforms.TaylorF2}
        @test s isa WaveformSettings{GWWaveforms.IMRPhenomD}
        @test_throws ArgumentError WaveformSettings("IMRPhenomXPHM")
        @test_throws ArgumentError WaveformSettings("IMRPhenomD"; mode_array = [(2, 2)])
        @test_throws ArgumentError WaveformSettings("TaylorF2"; pn_amplitude_order = 1)
        @test_throws ArgumentError WaveformSettings("TaylorF2"; pn_phase_order = 8)
        @test_throws ArgumentError WaveformSettings("TaylorF2"; pn_spin_order = -2)
        @test sprint(show, WaveformSettings("TaylorF2"; maximum_frequency = 300)) ==
            "WaveformSettings(\"TaylorF2\"; reference_frequency=50.0, minimum_frequency=20.0, maximum_frequency=300.0)"
    end

    @testset "parameters" begin
        s = WaveformSettings("IMRPhenomD")
        @test_throws ArgumentError hphc_fd(f, merge(base, (; foo = 1.0)), s)
        @test_throws ArgumentError hphc_fd(f, Base.structdiff(base, (; phase = 0.0)), s)
        # Precessing spins and eccentricity are outside both models.
        @test_throws ArgumentError hphc_fd(f, merge(base, (; a_1 = 0.5, tilt_1 = 0.7)), s)
        @test_throws ArgumentError hphc_fd(f, merge(base, (; eccentricity = 0.1)), s)
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

    @testset "grid" begin
        s = WaveformSettings("IMRPhenomD")
        @test_throws ArgumentError hphc_fd(range(1, 512; step = 0.25), base, s)
        @test_throws ArgumentError hphc_fd(range(0, 512; step = -0.25), base, s)
        @test_throws ArgumentError hphc_fd(range(0, 0; length = 1), base, s)
    end

    @testset "masking and f_max" begin
        for name in ("IMRPhenomD", "TaylorF2")
            h = hphc_fd(f, base, WaveformSettings(name; minimum_frequency = 30, maximum_frequency = 200))
            below = f .< 30
            above = f .> 200
            inside = .!(below .| above)
            @test all(iszero, h.plus[below]) && all(iszero, h.cross[below])
            @test all(iszero, h.plus[above]) && all(iszero, h.cross[above])
            @test all(!iszero, h.plus[inside][1:(end - 1)])
            # `maximum_frequency = nothing` means the last frequency of the grid.
            h_default = hphc_fd(f, base, WaveformSettings(name; minimum_frequency = 30))
            h_explicit = hphc_fd(
                f, base, WaveformSettings(name; minimum_frequency = 30, maximum_frequency = last(f))
            )
            @test h_default.plus == h_explicit.plus
            @test h_default.cross == h_explicit.cross
        end
    end

    @testset "orientation" begin
        # Face-on: hc = -i hp; edge-on: no cross polarisation, hp halved.
        for name in ("IMRPhenomD", "TaylorF2")
            s = WaveformSettings(name)
            face = hphc_fd(f, merge(base, (; theta_jn = 0.0)), s)
            @test face.cross ≈ -im .* face.plus
            edge = hphc_fd(f, merge(base, (; theta_jn = π / 2)), s)
            @test maximum(abs, edge.cross) < 1.0e-12 * maximum(abs, edge.plus)
            @test edge.plus ≈ face.plus ./ 2
        end
    end

    @testset "element type" begin
        s = WaveformSettings("TaylorF2")
        @test eltype(hphc_fd(f, base, s).plus) === ComplexF64
        @test eltype(hphc_fd(f, merge(base, (; mass_1 = 36)), s).plus) === ComplexF64
    end
end
