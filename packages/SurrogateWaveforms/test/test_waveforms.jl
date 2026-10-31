using TestItems

@testsnippet Setup_waveforms begin
    using FiniteDifferences
    using ForwardDiff
    using JET
    using LinearAlgebra
    using NPZ
    using Random
    using AbstractWaveforms
    using AbstractWaveformsTests
end

@testitem "waveform assembly kernels" setup = [Common, Setup_waveforms] begin
    @testset "quat_multiply" begin
        q = (0.5, -0.5, 0.5, 0.5)
        qc = (q[1], -q[2], -q[3], -q[4])
        @test all(SurrogateWaveforms.quat_multiply(q, qc) .≈ (1.0, 0.0, 0.0, 0.0))
        a = (0.9, 0.1, -0.3, 0.2)
        b = (0.2, 0.6, 0.1, -0.7)
        c = (-0.3, 0.2, 0.8, 0.4)
        @test all(
            SurrogateWaveforms.quat_multiply(SurrogateWaveforms.quat_multiply(a, b), c) .≈
                SurrogateWaveforms.quat_multiply(a, SurrogateWaveforms.quat_multiply(b, c))
        )
    end

    @testset "rotate_modes" begin
        rng = MersenneTwister(5)
        n = 4
        h = randn(rng, ComplexF64, n, SurrogateWaveforms.n_modes(4))
        # The identity rotor leaves the modes unchanged.
        ident = repeat([1.0 0.0 0.0 0.0], n)
        @test SurrogateWaveforms.rotate_modes(h, ident, 4) ≈ h
        # A rotation preserves the norm of each ell block.
        quat = randn(rng, n, 4)
        quat ./= sqrt.(sum(abs2, quat; dims = 2))
        r = SurrogateWaveforms.rotate_modes(h, quat, 4)
        for l in 2:4, k in 1:n
            cols = [SurrogateWaveforms.mode_column(l, m) for m in (-l):l]
            @test norm(r[k, cols]) ≈ norm(h[k, cols])
        end
        # Rotating by R then by its inverse recovers the modes.
        qinv = copy(quat)
        qinv[:, 2:4] .*= -1
        @test SurrogateWaveforms.rotate_modes(r, qinv, 4) ≈ h
    end

    @testset "rotate_modes against the full Wigner matrices" begin
        rng = MersenneTwister(6)
        n = 5
        function reference(h, quat, ell_max)
            out = zeros(promote_type(eltype(h), eltype(quat)), size(h))
            for k in 1:n, l in 2:ell_max
                D = SurrogateWaveforms.wigner_D(l, quat[k, 1], quat[k, 2], quat[k, 3], quat[k, 4])
                for mp in (-l):l, m in (-l):l
                    out[k, SurrogateWaveforms.mode_column(l, mp)] += D[mp + l + 1, m + l + 1] * h[k, SurrogateWaveforms.mode_column(l, m)]
                end
            end
            return out
        end
        for ell_max in 2:5
            h = randn(rng, ComplexF64, n, SurrogateWaveforms.n_modes(ell_max))
            quat = randn(rng, n, 4)
            quat ./= sqrt.(sum(abs2, quat; dims = 2))
            @test SurrogateWaveforms.rotate_modes(h, quat, ell_max) ≈ reference(h, quat, ell_max) rtol = 1.0e-13
        end
        # Dual numbers in the rotor and in the modes.
        h = randn(rng, ComplexF64, n, SurrogateWaveforms.n_modes(4))
        quat0 = randn(rng, n, 4)
        quat0 ./= sqrt.(sum(abs2, quat0; dims = 2))
        f(rotate) = x -> sum(abs2, rotate(h .* x[1], quat0 .+ x[2] .* [0.1 0.2 -0.1 0.05], 4))
        x0 = [1.1, 0.3]
        @test ForwardDiff.gradient(f(SurrogateWaveforms.rotate_modes), x0) ≈ ForwardDiff.gradient(f(reference), x0)
        @test_throws ArgumentError SurrogateWaveforms.rotate_modes(h, quat0, 6)
    end

    @testset "rotate_modes! reuses the buffer" begin
        rng = MersenneTwister(7)
        n = 6
        h = randn(rng, ComplexF64, n, SurrogateWaveforms.n_modes(5))
        quat = randn(rng, n, 4)
        quat ./= sqrt.(sum(abs2, quat; dims = 2))
        expected = SurrogateWaveforms.rotate_modes(h, quat, 5)
        buffer = copy(h)
        result = SurrogateWaveforms.rotate_modes!(buffer, quat, 5)
        @test result === buffer
        @test result == expected
        # A rotor of a wider element type cannot be stored in `h`: the result is a new array.
        wide = SurrogateWaveforms.rotate_modes!(copy(h), ForwardDiff.Dual.(quat, 1.0), 5)
        @test eltype(wide) <: Complex{<:ForwardDiff.Dual}
        @test ForwardDiff.value.(real.(wide)) ≈ real.(expected)
    end

    @testset "shared harmonics against retained Wigner-D" begin
        for theta in (0.0, 1.0e-10, 0.4, 1.7, π - 1.0e-10, Float64(π)), phi in (0.0, 0.7)
            for l in 2:5, m in (-l):l
                d = real(SurrogateWaveforms.wigner_D_element(l, m, 2, complex(cos(theta / 2)), complex(sin(theta / 2))))
                expected = sqrt((2l + 1) / (4π)) * d * cis(m * phi)
                @test SurrogateWaveforms.sylm_m2(l, m, theta, phi) ≈ expected atol = 3.0e-15 rtol = 2.0e-13
            end
        end
    end

    @testset "mode_sum on a single mode" begin
        modes = zeros(ComplexF64, 1, SurrogateWaveforms.n_modes(2))
        modes[1, SurrogateWaveforms.mode_column(2, 2)] = 1.0
        # ₋₂Y₂₂(θ, φ) = sqrt(5 / 64π) (1 + cos θ)² e^{2iφ}, with φ = π/2 - phi_ref.
        theta, phi_ref = 0.9, 0.4
        expected = sqrt(5 / (64π)) * (1 + cos(theta))^2 * cis(2 * (π / 2 - phi_ref))
        @test only(SurrogateWaveforms.mode_sum(modes, 2, theta, phi_ref)) ≈ expected
    end
end

@testitem "precessing waveform" tags = [:reference] setup = [Common, Setup_waveforms] begin
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
            sur = SurrogateWaveforms.load_surrogate(model)
            @test SurrogateWaveforms.ell_max(sur) == (model === :NRSur7dq4 ? 4 : 5)

            @testset "spin-weighted harmonics" begin
                for (i, (th, ph)) in enumerate(zip(ref["sylm_thetas"], ref["sylm_phis"]))
                    j = 0
                    for l in 2:5, m in (-l):l
                        j += 1
                        @test SurrogateWaveforms.sylm_m2(l, m, th, ph) ≈ ref["sylm"][i, j] atol = 1.0e-14
                    end
                end
            end

            @testset "modes, case $i" for i in 0:(ref["n_cases"] - 1)
                dt = ref["case$(i)_dt"]
                times = get(ref, "case$(i)_times", nothing)
                out = SurrogateWaveforms.inertial_modes(
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
                gen = SurrogateWaveforms.PrecessingGenerator(sur; f_ref = ref["mks_f_ref"])
                out = SurrogateWaveforms.hphc_td(gen, params; dt = ref["mks_dt"])
                grid = times(out)[2:(end - 1)]
                explicit = SurrogateWaveforms.hphc_td(grid, params, gen)
                prototype = TimeSeries(zeros(length(grid)); epoch = first(grid), sampling_interval = step(grid), channels = (:unused,))
                sampled = SurrogateWaveforms.hphc_td(prototype, params, gen)
                @test explicit isa Matrix
                @test explicit ≈ parent(out)[2:(end - 1), :]
                @test parent(sampled) == explicit
                @test times(sampled) == grid
                @test channels(sampled) == (:plus, :cross)
                irregular = collect(grid)[[1, 3, 8, 9]]
                @test SurrogateWaveforms.hphc_td(irregular, params, gen) ≈ explicit[[1, 3, 8, 9], :]
                modes_si = SurrogateWaveforms.waveform_modes_td(gen, params; dt = ref["mks_dt"])
                @test all(isapprox.(modes_si.t, times(out); rtol = 0, atol = 8eps(maximum(abs, modes_si.t))))
                @test modes_si.mode_indices == Tuple((l, m) for l in 2:gen.ell_max for m in (-l):l)
                unseen = merge(params, (; inclination = 0.1, phi_ref = 0.8))
                test_modes_td_contract(gen; params, viewing_params = unseen, dt = ref["mks_dt"])
                @test SurrogateWaveforms.waveform_modes_td(gen, unseen; dt = ref["mks_dt"]) == modes_si
                @test SurrogateWaveforms.AbstractWaveforms.approximant(gen) == String(model)
                @test SurrogateWaveforms.AbstractWaveforms.spin_support(gen) == :precessing
                @test SurrogateWaveforms.AbstractWaveforms.is_differentiable(gen)
                @test SurrogateWaveforms.AbstractWaveforms.supports_mode_array(gen)
                @test reference_frequency(gen) == gen.f_ref
                @test minimum_frequency(gen) == gen.f_low
                @test maximum_frequency(gen) === nothing

                h = out[Chan = At(:plus)] .- im .* out[Chan = At(:cross)]
                expected = ref["mks_h"]
                @test length(h) == length(expected)
                keep = collect(times(out)) .<= sur.coorbital.t[end - 1] * SurrogateWaveforms.MSUN_IN_SECONDS * ref["mks_M"]
                @test maximum(abs, collect(times(out)) - ref["mks_t"]) < 1.0e-12
                @test maximum(abs, h[keep] - expected[keep]) < 1.0e-11 * maximum(abs, expected)
            end

            @testset "dt grid equals the explicit grid" begin
                a = (0.3, 0.2, 0.4)
                b = (-0.1, 0.5, -0.3)
                by_dt = SurrogateWaveforms.inertial_modes(sur, 2.0, a, b; ell_max = 2, t_ref = -100.0, dt = 7.0)
                t0 = sur.coorbital.t[1]
                @test by_dt.t[1] == t0
                @test all(diff(by_dt.t) .≈ 7.0)
                @test by_dt.t[end] < sur.coorbital.t[end]
                @test by_dt.t[end] + 7.0 >= sur.coorbital.t[end]
                by_times = SurrogateWaveforms.inertial_modes(sur, 2.0, a, b; ell_max = 2, t_ref = -100.0, times = by_dt.t)
                @test by_times.modes == by_dt.modes
            end

            @testset "argument checks" begin
                a = (0.3, 0.2, 0.4)
                @test_throws ArgumentError SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; ell_max = 1)
                @test_throws ArgumentError SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; dt = 1.0, times = [0.0])
                @test_throws ArgumentError SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; t_ref = -100.0, f_ref = 0.01)
                @test_throws ArgumentError SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; times = [-6000.0, 0.0])
                @test_throws ArgumentError SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; times = [-100.0, 200.0])
                @test_throws ArgumentError SurrogateWaveforms.PrecessingGenerator(sur; f_ref = 10.0, f_low = 20.0)
            end

            @testset "ForwardDiff agrees with finite differences" begin
                fdm = central_fdm(5, 1)
                gen = SurrogateWaveforms.PrecessingGenerator(sur; ell_max = 3, f_ref = 25.0)
                function objective(p)
                    params = (
                        total_mass = p[1], luminosity_distance = 400.0, mass_ratio = p[2],
                        spin_1 = (p[3], p[4], p[5]), spin_2 = (p[6], 0.5, -0.3), inclination = p[7], phi_ref = p[8],
                    )
                    out = SurrogateWaveforms.hphc_td(gen, params; times = range(-0.4, 0.01; length = 40))
                    return 1.0e21 * sum(out[Chan = At(:plus)] .* range(1, 2; length = 40) .+ 2 .* out[Chan = At(:cross)])
                end
                p0 = [60.0, 2.2, 0.3, 0.2, 0.4, -0.1, 0.7, 0.3]
                g = ForwardDiff.gradient(objective, p0)
                @test all(isfinite, g)
                @test g ≈ grad(fdm, objective, p0)[1] rtol = 1.0e-5 atol = 1.0e-5
            end

            if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
                @testset "JET optimisation" begin
                    a = (0.3, 0.2, 0.4)
                    JET.@test_opt target_modules = (SurrogateWaveforms,) SurrogateWaveforms.inertial_modes(sur, 2.0, a, a)
                    JET.@test_opt target_modules = (SurrogateWaveforms,) SurrogateWaveforms.inertial_modes(sur, 2.0, a, a; f_ref = 0.015, dt = 2.0)
                end
            end
        end
    end
end
