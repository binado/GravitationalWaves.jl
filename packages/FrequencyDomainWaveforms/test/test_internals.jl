import ForwardDiff as FD

@testset "model building blocks" begin
    G = FrequencyDomainWaveforms

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
        for Mf in (0.005, 0.02, 0.05, 0.1)
            @test isfinite(G.phenomd_amplitude(m, Mf))
            @test isfinite(G.phenomd_phase(m, Mf))
        end

        # The amplitude is continuous at the region boundaries, the phase and its derivative too.
        ε = 1.0e-9
        for Mf in (G.AMP_FJOIN_INS, m.amp.f_max_calc)
            @test G.phenomd_amplitude(m, Mf - ε) ≈ G.phenomd_amplitude(m, Mf + ε) rtol = 1.0e-6
        end
        for Mf in (G.PHI_FJOIN_INS, m.conn.f_mrd_join)
            @test G.phenomd_phase(m, Mf - ε) ≈ G.phenomd_phase(m, Mf + ε) rtol = 1.0e-6
            d_left = (G.phenomd_phase(m, Mf - ε) - G.phenomd_phase(m, Mf - 2ε)) / ε
            d_right = (G.phenomd_phase(m, Mf + 2ε) - G.phenomd_phase(m, Mf + ε)) / ε
            @test d_left ≈ d_right rtol = 1.0e-3
        end
        # The Lorentzian-like merger-ringdown ansatz peaks at fmax_calc.
        peak = G.amp_mrd_ansatz(m.amp.f_max_calc, m.amp)
        @test peak > G.amp_mrd_ansatz(0.99 * m.amp.f_max_calc, m.amp)
        @test peak > G.amp_mrd_ansatz(1.01 * m.amp.f_max_calc, m.amp)

        # Eta is snapped only by rounding.
        @test G.nudged_eta(1.0, 1.0) == 0.25
        @test G.nudged_eta(3.0, 1.0) == 3.0 / 16.0
        @test G.nudged_eta(1.0, 1.0 + 1.0e-16) == 0.25
        # `fmax_calc` takes the second branch for gamma2 > 1 without evaluating a NaN square root.
        @test G.fmax_calc(0.07, 0.01, 1.5, 1.0) == abs(0.07 + (-0.01 * 1.0) / 1.5)
        @test G.fmax_calc(0.07, 0.01, 0.5, 1.0) == abs(0.07 + (0.01 * (-1 + sqrt(1 - 0.25)) * 1.0) / 0.5)

        # Equal masses sit at the upper limit of eta.
        @test isfinite(G.phenomd_amplitude(G.phenomd_model(30.0, 30.0, 0.0, 0.0, 50.0, 0.0), 0.01))
    end

    @testset "all regions stay finite where they are not selected" begin
        # The kernels evaluate the inspiral, intermediate and merger-ringdown pieces everywhere and
        # select one, so each must be finite (value and derivative) for every Mf > 0 it may see.
        for (m1, m2) in ((36.0, 29.0), (10.0, 1.4), (200.0, 150.0), (3000.0, 3000.0), (50.0, 1.0)),
                (chi1, chi2) in ((0.0, 0.0), (0.99, 0.99), (-0.99, -0.99), (0.9, -0.7), (-0.5, 0.95))

            m = G.phenomd_model(m1, m2, chi1, chi2, 20.0, 0.3)
            pieces = (
                Mf -> G.amp_ins_ansatz(Mf, G.useful_powers(Mf), m.amp.pre),
                Mf -> G.amp_int_ansatz(Mf, m.amp),
                Mf -> G.amp_mrd_ansatz(Mf, m.amp),
                Mf -> G.phi_ins_ansatz_int(Mf, G.useful_powers(Mf), m.phi_pre, m.phi),
                Mf -> m.phi.eta_inv * G.phi_int_ansatz(Mf, m.phi) + m.conn.c1_int + m.conn.c2_int * Mf,
                Mf -> m.phi.eta_inv * G.phi_mrd_ansatz_int(Mf, m.phi) + m.conn.c1_mrd + m.conn.c2_mrd * Mf,
            )
            for piece in pieces, Mf in exp10.(range(-4, log10(10.0); length = 60))
                @test isfinite(piece(Mf))
                @test isfinite(FD.derivative(piece, Mf))
            end
        end
    end

    @testset "gradients across the region boundaries" begin
        m = G.phenomd_model(36.0, 29.0, 0.5, 0.3, 50.0, 1.3)
        ε = 1.0e-7
        for edge in (G.AMP_FJOIN_INS, m.amp.f_max_calc), Mf in (edge - ε, edge, edge + ε)
            @test isfinite(FD.derivative(x -> G.phenomd_amplitude(m, x), Mf))
        end
        # The phase is C¹ across both joins.
        for edge in (G.PHI_FJOIN_INS, m.conn.f_mrd_join)
            left = FD.derivative(x -> G.phenomd_phase(m, x), edge - ε)
            right = FD.derivative(x -> G.phenomd_phase(m, x), edge + ε)
            @test isfinite(left) && isfinite(right)
            @test left ≈ right rtol = 1.0e-3
        end
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

    @testset "spline evaluation is bit-identical to the searchsorted form" begin
        # The previous implementation: binary search, then the coefficients from the raw arrays.
        # A natural spline has c[n] = 0, so the full coefficient vector is recovered from the segments.
        function reference_eval(x_knots, y_knots, c, x)
            n = length(x_knots)
            i = clamp(searchsortedlast(x_knots, x), 1, n - 1)
            x_lo = x_knots[i]
            dx = x_knots[i + 1] - x_lo
            y_lo = y_knots[i]
            dy = y_knots[i + 1] - y_lo
            delx = x - x_lo
            c_i = c[i]
            c_ip1 = c[i + 1]
            b_i = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0
            d_i = (c_ip1 - c_i) / (3.0 * dx)
            return y_lo + delx * (b_i + delx * (c_i + delx * d_i))
        end
        for (spline, y) in ((G.FRING_SPLINE, G.QNM_FRING), (G.FDAMP_SPLINE, G.QNM_FDAMP))
            knots = spline.x
            c = vcat(spline.c, 0.0)
            mids = (knots[1:(end - 1)] .+ knots[2:end]) ./ 2
            rng = range(-1.2, 1.2; length = 97)
            for x in Iterators.flatten((knots, mids, rng, (-1.0, 1.0)))
                @test G.spline_eval(spline, x) == reference_eval(knots, y, c, x)
            end
        end
    end
end

@testset "TaylorF2 amplitude_phase_fd is zero outside the band" begin
    f = range(0, 256; step = 0.25)
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
        a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
    )
    for g in (
            TaylorF2(; minimum_frequency = 25.0),
            TaylorF2(; minimum_frequency = 25.0, maximum_frequency = 180.0, pn_amplitude_order = 3),
            TaylorF2(; minimum_frequency = 20.0, pn_phase_order = 4, pn_spin_order = 3, reference_frequency = 0),
        )
        f_min, f_max = FrequencyDomainWaveforms.frequency_bounds(g, f)
        ap = amplitude_phase_fd(f, params, g)
        inband = [f_min <= fi <= f_max for fi in f]
        @test all(iszero, ap.amplitude[.!inband]) && all(iszero, ap.phase[.!inband])
        @test all(>(0), ap.amplitude[inband])
    end
    # Nothing in band.
    @test all(iszero, amplitude_phase_fd(f, params, TaylorF2(; minimum_frequency = 300.0, maximum_frequency = 400.0)).amplitude)
end

@testset "IMRPhenomD requires mass_1 ≥ mass_2" begin
    f = range(0, 512; step = 0.25)
    g = IMRPhenomD(; minimum_frequency = 20.0, maximum_frequency = 500.0)
    base = (; luminosity_distance = 400.0, theta_jn = 0.4, phase = 1.2)
    a = hphc_fd(f, merge(base, (; mass_1 = 36.0, mass_2 = 29.0, a_1 = 0.5, a_2 = 0.3)), g)
    @test isfinite(sum(abs2, a[:, :plus]))
    @test_throws ArgumentError validate_parameters(
        f, merge(base, (; mass_1 = 29.0, mass_2 = 36.0, a_1 = 0.3, a_2 = 0.5)), g,
    )
end

@testset "IMRPhenomD amplitude_phase_fd band and cutoff" begin
    f = range(0, 256; step = 0.25)
    base = (
        luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
        a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
    )
    df = Float64(step(f))
    fcut(m1, m2) = FrequencyDomainWaveforms.F_CUT / ((m1 + m2) * FrequencyDomainWaveforms.MTSUN_SI)
    # The last pair puts f_cut inside the grid; the second-last makes it an exact grid point.
    for (m1, m2) in ((36.0, 29.0), (29.0, 36.0), (90.0, 85.0), (130.0, 100.0))
        for g in (
                IMRPhenomD(; minimum_frequency = 25.0),
                IMRPhenomD(; minimum_frequency = 25.0, maximum_frequency = 180.0),
                IMRPhenomD(; minimum_frequency = 20.1, maximum_frequency = 200.3, reference_frequency = 0),
                IMRPhenomD(; minimum_frequency = 25.0, maximum_frequency = fcut(m1, m2)),
            )
            f_min, f_max = FrequencyDomainWaveforms.frequency_bounds(g, f)
            fcut(m1, m2) > f_min || continue
            ap = amplitude_phase_fd(f, merge(base, (; mass_1 = m1, mass_2 = m2)), g)
            # LAL excludes the rounded-down upper edge.
            upper = floor(min(fcut(m1, m2), f_max) / df) * df
            for (i, fi) in enumerate(f)
                if f_min <= fi < upper
                    @test ap.amplitude[i] > 0
                else
                    @test iszero(ap.amplitude[i]) && iszero(ap.phase[i])
                end
            end
        end
    end
    # Nothing in band.
    g = IMRPhenomD(; minimum_frequency = 300.0, maximum_frequency = 400.0)
    @test all(iszero, amplitude_phase_fd(f, merge(base, (; mass_1 = 36.0, mass_2 = 29.0)), g).amplitude)
end
