using FiniteDifferences
using ForwardDiff
using JET
using LinearAlgebra
using NPZ
using Random

@testset "waveform assembly kernels" begin
    @testset "quat_multiply" begin
        q = (0.5, -0.5, 0.5, 0.5)
        qc = (q[1], -q[2], -q[3], -q[4])
        @test all(GWSurrogates.quat_multiply(q, qc) .≈ (1.0, 0.0, 0.0, 0.0))
        a = (0.9, 0.1, -0.3, 0.2)
        b = (0.2, 0.6, 0.1, -0.7)
        c = (-0.3, 0.2, 0.8, 0.4)
        @test all(
            GWSurrogates.quat_multiply(GWSurrogates.quat_multiply(a, b), c) .≈
                GWSurrogates.quat_multiply(a, GWSurrogates.quat_multiply(b, c))
        )
    end

    @testset "rotate_modes" begin
        rng = MersenneTwister(5)
        n = 4
        h = randn(rng, ComplexF64, n, GWSurrogates.n_modes(4))
        # The identity rotor leaves the modes unchanged.
        ident = repeat([1.0 0.0 0.0 0.0], n)
        @test GWSurrogates.rotate_modes(h, ident, 4) ≈ h
        # A rotation preserves the norm of each ell block.
        quat = randn(rng, n, 4)
        quat ./= sqrt.(sum(abs2, quat; dims = 2))
        r = GWSurrogates.rotate_modes(h, quat, 4)
        for l in 2:4, k in 1:n
            cols = [GWSurrogates.mode_column(l, m) for m in (-l):l]
            @test norm(r[k, cols]) ≈ norm(h[k, cols])
        end
        # Rotating by R then by its inverse recovers the modes.
        qinv = copy(quat)
        qinv[:, 2:4] .*= -1
        @test GWSurrogates.rotate_modes(r, qinv, 4) ≈ h
    end

    @testset "mode_sum on a single mode" begin
        modes = zeros(ComplexF64, 1, GWSurrogates.n_modes(2))
        modes[1, GWSurrogates.mode_column(2, 2)] = 1.0
        # ₋₂Y₂₂(θ, φ) = sqrt(5 / 64π) (1 + cos θ)² e^{2iφ}, with φ = π/2 - phi_ref.
        theta, phi_ref = 0.9, 0.4
        expected = sqrt(5 / (64π)) * (1 + cos(theta))^2 * cis(2 * (π / 2 - phi_ref))
        @test only(GWSurrogates.mode_sum(modes, 2, theta, phi_ref)) ≈ expected
    end
end

@testset "precessing waveform" begin
    ci = get(ENV, "CI", "false") == "true"
    for (model, params) in pairs(MODEL_PARAMS)
        fixture = joinpath(FIXTURE_DIR, "$(model)_waveforms.npz")
        if !isfile(fixture)
            msg = "waveform reference fixture $fixture not found; generate it with python/generate_surrogates.py"
            ci ? error(msg) : @warn msg
            continue
        end
        @testset "$model" begin
            ref = npzread(fixture)
            sur = GWSurrogates.load_surrogate(model)
            @test GWSurrogates.ell_max(sur) == (model === :NRSur7dq4 ? 4 : 5)

            @testset "spin-weighted harmonics" begin
                for (i, (th, ph)) in enumerate(zip(ref["sylm_thetas"], ref["sylm_phis"]))
                    j = 0
                    for l in 2:5, m in (-l):l
                        j += 1
                        @test GWSurrogates.sylm_m2(l, m, th, ph) ≈ ref["sylm"][i, j] atol = 1.0e-14
                    end
                end
            end

            @testset "modes, case $i" for i in 0:(ref["n_cases"] - 1)
                dt = ref["case$(i)_dt"]
                times = get(ref, "case$(i)_times", nothing)
                out = GWSurrogates.inertial_modes(
                    sur, ref["case$(i)_q"], ref["case$(i)_chiA0"], ref["case$(i)_chiB0"];
                    ell_max = ref["case$(i)_ell_max"], f_ref = ref["case$(i)_f_ref"], f_low = ref["case$(i)_f_low"],
                    dt = isnan(dt) ? nothing : dt, times = times,
                    init_orbphase = ref["case$(i)_orbphase0"], init_quat = ref["case$(i)_quat0"],
                )
                expected = permutedims(ref["case$(i)_modes"])
                @test size(out.modes) == size(expected)
                @test maximum(abs, out.t - ref["case$(i)_t"]) < 1.0e-9
                # gwsurrogate's C++ spline mis-evaluates isolated points in the last interval of
                # the model grid; this package evaluates them correctly, so compare elsewhere.
                keep = out.t .<= sur.coorbital.t[end - 1]
                @test count(keep) >= length(keep) - 2
                @test maximum(abs, out.modes[keep, :] - expected[keep, :]) < 1.0e-11 * maximum(abs, expected)
            end

            @testset "SI strain end to end" begin
                params = (
                    total_mass = ref["mks_M"], luminosity_distance = ref["mks_dist_mpc"], mass_ratio = ref["mks_q"],
                    spin_1 = ref["mks_chiA"], spin_2 = ref["mks_chiB"], inclination = ref["mks_inclination"],
                    phi_ref = ref["mks_phi_ref"],
                )
                gen = GWSurrogates.PrecessingGenerator(sur; f_ref = ref["mks_f_ref"])
                out = GWSurrogates.hphc_td(gen, params; dt = ref["mks_dt"])
                h = out.hplus .- im .* out.hcross
                expected = ref["mks_h"]
                @test length(h) == length(expected)
                keep = out.t .<= sur.coorbital.t[end - 1] * GWSurrogates.MSUN_IN_SECONDS * ref["mks_M"]
                @test maximum(abs, out.t - ref["mks_t"]) < 1.0e-12
                @test maximum(abs, h[keep] - expected[keep]) < 1.0e-11 * maximum(abs, expected)
            end

            @testset "argument checks" begin
                a = (0.3, 0.2, 0.4)
                @test_throws ArgumentError GWSurrogates.inertial_modes(sur, 2.0, a, a; ell_max = 1)
                @test_throws ArgumentError GWSurrogates.inertial_modes(sur, 2.0, a, a; dt = 1.0, times = [0.0])
                @test_throws ArgumentError GWSurrogates.inertial_modes(sur, 2.0, a, a; t_ref = -100.0, f_ref = 0.01)
                @test_throws ArgumentError GWSurrogates.inertial_modes(sur, 2.0, a, a; times = [-6000.0, 0.0])
                @test_throws ArgumentError GWSurrogates.inertial_modes(sur, 2.0, a, a; times = [-100.0, 200.0])
                @test_throws ArgumentError GWSurrogates.PrecessingGenerator(sur; f_ref = 10.0, f_low = 20.0)
            end

            @testset "ForwardDiff agrees with finite differences" begin
                fdm = central_fdm(5, 1)
                gen = GWSurrogates.PrecessingGenerator(sur; ell_max = 3, f_ref = 25.0)
                function objective(p)
                    params = (
                        total_mass = p[1], luminosity_distance = 400.0, mass_ratio = p[2],
                        spin_1 = (p[3], p[4], p[5]), spin_2 = (p[6], 0.5, -0.3), inclination = p[7], phi_ref = p[8],
                    )
                    out = GWSurrogates.hphc_td(gen, params; times = range(-0.4, 0.01; length = 40))
                    return 1.0e21 * sum(out.hplus .* range(1, 2; length = 40) .+ 2 .* out.hcross)
                end
                p0 = [60.0, 2.2, 0.3, 0.2, 0.4, -0.1, 0.7, 0.3]
                g = ForwardDiff.gradient(objective, p0)
                @test all(isfinite, g)
                @test g ≈ grad(fdm, objective, p0)[1] rtol = 1.0e-5 atol = 1.0e-5
            end

            if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
                @testset "JET optimisation" begin
                    a = (0.3, 0.2, 0.4)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.inertial_modes(sur, 2.0, a, a)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.inertial_modes(sur, 2.0, a, a; f_ref = 0.015, dt = 2.0)
                end
            end
        end
    end
end
