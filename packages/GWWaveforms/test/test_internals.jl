@testset "model building blocks" begin
    G = GWWaveforms

    @testset "TaylorF2 phasing" begin
        eta = 36.0 * 29.0 / 65.0^2
        pn = G.taylorf2_pn_phasing(36.0, 29.0, 0.0, 0.0)
        @test pn.v[1] ≈ 3 / (128 * eta)
        @test pn.v[2] == 0
        # Without spin terms the spins do not matter.
        @test G.taylorf2_pn_phasing(36.0, 29.0, 0.5, -0.3; spin_order = 0) == pn
        @test G.taylorf2_pn_phasing(36.0, 29.0, 0.5, -0.3) != pn
        # Truncation to the Newtonian order keeps only the leading term.
        newtonian = G.truncate_phasing(G.taylorf2_pn_phasing(36.0, 29.0, 0.5, 0.3), 0)
        @test G.taylorf2_phasing(newtonian, 0.2) ≈ newtonian.v[1] / 0.2^5
        # The higher orders make a small, positive correction to the Newtonian phasing.
        full = G.taylorf2_pn_phasing(36.0, 29.0, 0.0, 0.0)
        @test 0 < abs(G.taylorf2_phasing(full, 0.18) / (full.v[1] / 0.18^5) - 1) < 0.2
    end

    @testset "IMRPhenomD" begin
        m = G.phenomd_model(36.0, 29.0, 0.5, 0.3, 50.0, 1.3)
        swapped = G.phenomd_model(29.0, 36.0, 0.3, 0.5, 50.0, 1.3)
        for Mf in (0.005, 0.02, 0.05, 0.1)
            @test G.phenomd_amplitude(m, Mf) == G.phenomd_amplitude(swapped, Mf)
            @test G.phenomd_phase(m, Mf) == G.phenomd_phase(swapped, Mf)
        end

        # The amplitude is continuous at the region boundaries, the phase and its derivative too.
        ε = 1.0e-9
        for Mf in (G.AMP_FJOIN_INS, m.amp.fmaxCalc)
            @test G.phenomd_amplitude(m, Mf - ε) ≈ G.phenomd_amplitude(m, Mf + ε) rtol = 1.0e-6
        end
        for Mf in (G.PHI_FJOIN_INS, m.conn.fMRDJoin)
            @test G.phenomd_phase(m, Mf - ε) ≈ G.phenomd_phase(m, Mf + ε) rtol = 1.0e-6
            d_left = (G.phenomd_phase(m, Mf - ε) - G.phenomd_phase(m, Mf - 2ε)) / ε
            d_right = (G.phenomd_phase(m, Mf + 2ε) - G.phenomd_phase(m, Mf + ε)) / ε
            @test d_left ≈ d_right rtol = 1.0e-3
        end
        # The Lorentzian-like merger-ringdown ansatz peaks at fmaxCalc.
        peak = G.amp_mrd_ansatz(m.amp.fmaxCalc, m.amp)
        @test peak > G.amp_mrd_ansatz(0.99 * m.amp.fmaxCalc, m.amp)
        @test peak > G.amp_mrd_ansatz(1.01 * m.amp.fmaxCalc, m.amp)

        # Equal masses sit at the upper limit of eta.
        @test isfinite(G.phenomd_amplitude(G.phenomd_model(30.0, 30.0, 0.0, 0.0, 50.0, 0.0), 0.01))
    end

    @testset "natural cubic spline" begin
        x = collect(range(0, 1; length = 41))
        s = G.NaturalCubicSpline(x, sin.(2 .* x))
        @test G.spline_eval(s, 0.3141) ≈ sin(2 * 0.3141) atol = 1.0e-5
        @test G.spline_eval(s, 0.0) == 0.0
        @test G.spline_eval(s, 1.0) ≈ sin(2.0) atol = 1.0e-12
        # The QNM tables: M ω = 0.37367 for the (2,2) mode of a Schwarzschild black hole.
        @test G.spline_eval(G.FRING_SPLINE, 0.0) ≈ 0.37367 / 2π atol = 1.0e-6
    end
end

@testset "TaylorF2 static slice against the loop" begin
    f = range(0, 256; step = 0.25)
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
        a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
    )
    p = GWWaveforms.aligned_parameters(params)
    df = Float64(step(f))
    for g in (
            TaylorF2(; minimum_frequency = 25.0),
            TaylorF2(; minimum_frequency = 25.0, maximum_frequency = 180.0, pn_amplitude_order = 3),
            TaylorF2(; minimum_frequency = 20.0, pn_phase_order = 4, pn_spin_order = 3, reference_frequency = 0),
        )
        f_min, f_max = GWWaveforms.frequency_bounds(g, f)
        f_start = g.pn_amplitude_order != 0 ?
            GWWaveforms.f_low_to_f_start(f_min, g.pn_amplitude_order) : f_min
        lo, hi = GWWaveforms.taylorf2_slice(f, f_start, f_min, f_max)
        expected = GWWaveforms.taylorf2_fd(Float64, length(f), df, p, g, f_start, f_max)
        expected[[!(f_min <= fi <= f_max) for fi in f]] .= 0
        got = GWWaveforms.taylorf2_fd_slice(Float64, length(f), lo, hi, df, p, g)
        @test got ≈ expected rtol = 1.0e-12
        @test all(iszero, got[1:lo]) && all(iszero, got[(hi + 2):end])
    end
    # Nothing in band.
    @test all(iszero, GWWaveforms.taylorf2_fd_slice(Float64, length(f), GWWaveforms.taylorf2_slice(f, 300.0, 300.0, 400.0)..., df, p, TaylorF2()))
end
