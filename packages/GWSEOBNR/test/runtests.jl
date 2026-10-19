using Test, GWSEOBNR, GWParameterConversion, ForwardDiff, LinearAlgebra, TOML
const S = GWSEOBNR

@testset "configuration and envelope" begin
    g = SEOBNRv5HM()
    @test approximant(g) == "SEOBNRv5HM"
    @test minimum_frequency(g) == reference_frequency(g) == 20
    @test maximum_frequency(g) === nothing
    @test spin_support(g) == :aligned
    @test !supports_mode_array(g)
    @test !is_differentiable(g)
    f = 0:0.0625:2048
    p = (mass_1 = 40.0, mass_2 = 20.0, luminosity_distance = 500.0, theta_jn = 0.7, phase = 0.3)
    @test validate_parameters(f, p, g) === nothing
    @test_throws ArgumentError validate_parameters(range(0.0; step = Inf, length = 2), p, g)
    @test validate_parameters(f, merge(p, (; a_1 = 0.4, tilt_1 = π)), g) === nothing
    @test validate_parameters(f, merge(p, (; chirp_mass = GWParameterConversion.chirp_mass(40.0, 20.0))), g) === nothing
    for nt in (
            (; mass_1 = 10.0), (; mass_2 = 0.0), (; a_1 = 0.7), (; a_1 = 0.1, tilt_1 = 0.5),
            (; eccentricity = 0.01), (; luminosity_distance = 0.0), (; phase = Inf), (; bad_key = 1),
        )
        @test_throws ArgumentError validate_parameters(f, merge(p, nt), g)
    end
    @test_throws ArgumentError SEOBNRv5HM(reference_frequency = 30)
    @test_throws ArgumentError SEOBNRv5HM(reltol = 0)
    @test_throws ArgumentError SEOBNRv5HM(maximum_frequency = 10)
end

@testset "pinned pySEOBNR component fidelity" begin
    oracle = TOML.parsefile(joinpath(@__DIR__, "oracle.toml"))
    @test oracle["pyseobnr_version"] == "0.3.7"
    for case in oracle["cases"]
        c = S.hamiltonian_coefficients(case["q"], 1.0, case["chi1"], case["chi2"])
        q, p = (12.0, 0.3), (-0.005, 4.0)
        d = S.hamiltonian_dynamics(q, p, c)
        @test collect(S.hamiltonian(q, p, c)) ≈ case["hamiltonian"] rtol = 1.0e-9 atol = 1.0e-14
        @test collect(S.hamiltonian_gradient(q, p, c)) ≈ case["gradient"] rtol = 1.0e-9 atol = 1.0e-14
        @test vec(S.hamiltonian_hessian(q, p, c)) ≈ case["hessian"] rtol = 1.0e-9 atol = 1.0e-14
        @test collect(d) ≈ case["dynamics"] rtol = 1.0e-9 atol = 1.0e-14
        @test c.a6 ≈ case["a6"] rtol = 1.0e-9
        @test c.dSO ≈ case["dSO"] rtol = 1.0e-9
        ap, am = c.X_1 * c.chi_1 + c.X_2 * c.chi_2, c.X_1 * c.chi_1 - c.X_2 * c.chi_2
        @test S.NR_deltaT_NS(c.nu) + S.NR_deltaT(c.nu, ap, am) ≈ case["delta_t"] rtol = 1.0e-9
        omega_circ = S.orbital_frequency(q, (0.0, 4.0), c)
        flux = S.radiation_flux(4.0, d[4], omega_circ, c.nu * d[5], c, S.flux_coefficients(c))
        @test flux ≈ case["flux"] rtol = 1.0e-9
        fake = (; sol = t -> [q..., p...], c, coeffs = S.flux_coefficients(c), phi_ref = 0.0)
        @test S.inspiral_mode(fake, 0.0) ≈ complex(case["mode_real"], case["mode_imag"]) rtol = 1.0e-9
        ic = S.initial_conditions(c.X_1, c.X_2, c.chi_1, c.chi_2, case["initial_omega"])
        @test [ic[1], ic[4], ic[3]] ≈ case["initial_conditions"] rtol = 1.0e-9 atol = 1.0e-12
        @test S.bbh_final_mass_non_precessing_UIB2016(c.X_1, c.X_2, c.chi_1, c.chi_2) ≈ case["final_mass"] rtol = 1.0e-9
        @test S.bbh_final_spin_non_precessing_HBR2016(c.X_1, c.X_2, c.chi_1, c.chi_2, "M3J4") ≈ case["final_spin"] rtol = 1.0e-9
        for name in ("iv_habs", "iv_hdot", "iv_hdotdot", "iv_omega", "iv_omegadot", "mr_c1f", "mr_c2f", "mr_d1f", "mr_d2f")
            @test getfield(S, Symbol(name))(c) ≈ case[name] rtol = 1.0e-9 atol = 1.0e-14
        end
        # Check the analytic gradients/Hessian and generic scalar arithmetic independently.
        state = [q..., p...]
        @test ForwardDiff.gradient(x -> S.hamiltonian((x[1], x[2]), (x[3], x[4]), c)[1], state) ≈ collect(S.hamiltonian_gradient(q, p, c)) rtol = 1.0e-9 atol = 1.0e-13
        @test ForwardDiff.hessian(x -> S.hamiltonian((x[1], x[2]), (x[3], x[4]), c)[1], state) ≈ S.hamiltonian_hessian(q, p, c) rtol = 1.0e-9 atol = 1.0e-13
    end
end

@testset "bundled QNM interpolants" begin
    for case in TOML.parsefile(joinpath(@__DIR__, "qnm_oracle.toml"))["cases"]
        spin = case["m"] > 0 ? case["spin"] : -case["spin"]
        @test S.qnm_220(spin, 1.0) ≈ complex(case["real"], case["imag"]) rtol = 1.0e-10
    end
    @test_throws ArgumentError S.qnm_220(1.1, 1.0)
    @test S.qnm_220(0.7, 0.95) ≈ S.qnm_220(0.7, 1.0) / 0.95
end

@testset "local tail and FFT ForwardDiff adapters" begin
    for N in (8, 10)
        x = [sin(0.7j) + 0.2cos(0.3j) for j in 1:N]
        direct(x) = [sum(x[j + 1] * cis(-2π * k * j / N) for j in 0:(N - 1)) for k in 0:div(N, 2)]
        @test S.real_fft(x) ≈ direct(x) rtol = 1.0e-13 atol = 1.0e-13
        v = reverse(x)
        ad = ForwardDiff.derivative(t -> S.real_fft(x + t * v), 0.0)
        @test ad ≈ direct(v) rtol = 1.0e-13 atol = 1.0e-13
        dt, tstart = 0.2, -0.7
        @test S.physical_spectrum(x, dt, tstart) ≈ [dt * sum(x[j + 1] * cis(-2π * k * (tstart + j * dt) / (N * dt)) for j in 0:(N - 1)) for k in 0:div(N, 2)] rtol = 1.0e-13 atol = 1.0e-13
    end
    tail(t) = S.complex_tail(2, 2, 0.03 + 0.01t, 0.98 + 0.02t)
    eps = 1.0e-4
    @test ForwardDiff.derivative(tail, 0.0) ≈ (tail(eps) - tail(-eps)) / (2eps) rtol = 1.0e-8
    @test S.time_jet(t -> 2 + 3t + 4t^2, 0.7) ≈ [6.06, 8.6, 8.0] rtol = 1.0e-10
end

@testset "roots and explicit failures" begin
    @test S.check_evolution_retcode(S.SciMLBase.ReturnCode.Terminated) === nothing
    @test_throws "EOB evolution exhausted time coverage without a physical guard" S.check_evolution_retcode(S.SciMLBase.ReturnCode.Success)
    @test_throws "EOB evolution failed:" S.check_evolution_retcode(S.SciMLBase.ReturnCode.MaxIters)
    @test_throws S.SolveError S.check_isco_interval(10.0, 4.0, 3.0)
    @test_throws S.SolveError S.check_isco_interval(2.0, 1.5, 3.0)
    @test_throws S.SolveError S.attachment_time(5.0, -6.0, 0.0, 10.0)
    @test_throws S.SolveError S.attachment_time(5.0, 6.0, 0.0, 10.0)
    @test S.attachment_time(5.0, 1.0, 0.0, 10.0) == 4.0
    @test_throws S.SolveError S.bracket_root(x -> x^2 + 1, -1, 1)
    @test_throws S.SolveError S.continuous_root(x -> x^3, -1, 1, x -> 3x^2)
    @test_throws S.SolveError S.newton_root((x, p) -> [1.0], [1.0], [1.0])
    root(t) = only(S.implicit_root((x, p) -> [x[1]^2 - p[1]], [1.0], [t]))
    @test root(2.0) ≈ sqrt(2)
    @test ForwardDiff.derivative(root, 2.0) ≈ 1 / (2sqrt(2))
    timing(t) = S.continuous_root(s -> s^2 - t, 0.1, 3.0, s -> 2s)
    @test ForwardDiff.derivative(timing, 2.0) ≈ 1 / (2sqrt(2)) rtol = 1.0e-10
end

@testset "short native IMR and boundaries" begin
    f = 0:0.25:1024
    p = (mass_1 = 50.0, mass_2 = 50.0, luminosity_distance = 500.0, theta_jn = 0.7, phase = 0.3)
    g = SEOBNRv5HM(maximum_frequency = 800)
    h = hphc_fd(f, p, g)
    @test length(h.plus) == length(h.cross) == length(f)
    @test all(isfinite, h.plus) && all(isfinite, h.cross)
    @test all(iszero, h.plus[(f .< 20) .| (f .> 800)])
    @test maximum(abs, h.plus) > 0
    @test hphc_fd(f, merge(p, (; luminosity_distance = 1000.0)), g).plus ≈ h.plus / 2 rtol = 1.0e-12
    @test_throws S.SolveError hphc_fd(0:1.0:1024, p, g)
    @test_throws S.SolveError hphc_fd(f, p, SEOBNRv5HM(post_merger_duration = 0.1))
    sym(t) = sum(abs2, hphc_fd(f, merge(p, (; mass_1 = 50 - t, mass_2 = 50 - t)), g).plus)
    @test isfinite(ForwardDiff.derivative(sym, 0.0))
end

@testset "compact oracle waveform anchors" begin
    ref = TOML.parsefile(joinpath(@__DIR__, "short_waveform.toml"))
    f = 0:ref["df"]:ref["fmax"]
    p = (
        mass_1 = 32.0, mass_2 = 8.0, a_1 = 0.6, a_2 = 0.6, tilt_1 = π,
        luminosity_distance = 500.0, theta_jn = 0.7, phase = 0.3,
    )
    h = hphc_fd(f, p, SEOBNRv5HM())
    i = ref["indices"]
    expected_plus = complex.(ref["plus_real"], ref["plus_imag"])
    expected_cross = complex.(ref["cross_real"], ref["cross_imag"])
    @test norm(h.plus[i] - expected_plus) / norm(expected_plus) < 1.0e-3
    @test norm(h.cross[i] - expected_cross) / norm(expected_cross) < 1.0e-3
    @test norm(h.plus) ≈ ref["plus_norm"] rtol = 1.0e-3
    @test norm(h.cross) ≈ ref["cross_norm"] rtol = 1.0e-3
end
