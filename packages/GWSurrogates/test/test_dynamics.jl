using FiniteDifferences
using ForwardDiff
using JET
using NPZ
using StaticArrays

@testset "dynamics kernels" begin
    @testset "ab4_dy integrates a cubic derivative exactly" begin
        # Derivative samples k(t) = 1 + 2t - t^2 + 0.5t^3 at t = 0, dt1, dt1+dt2, dt1+dt2+dt3.
        k(t) = 1 + 2t - t^2 + 0.5t^3
        K(t) = t + t^2 - t^3 / 3 + t^4 / 8
        dt1, dt2, dt3, dt4 = 0.7, 1.1, 0.9, 1.3
        t4 = dt1 + dt2 + dt3
        ks = (k(0.0), k(dt1), k(dt1 + dt2), k(t4))
        incr = GWSurrogates.ab4_dy(ks..., dt1, dt2, dt3, dt4)
        @test incr ≈ K(t4 + dt4) - K(t4) rtol = 1.0e-12
        # It acts component-wise on state vectors.
        v = GWSurrogates.ab4_dy(map(x -> SVector(x, 2x), ks)..., dt1, dt2, dt3, dt4)
        @test v ≈ SVector(incr, 2incr)
    end

    @testset "interp4 is a natural cubic spline" begin
        ts = (93.8, 95.8, 97.9, 100.0)
        us = (SVector(0.1, 1.0), SVector(0.5, -1.0), SVector(0.2, 0.3), SVector(0.9, 2.0))
        sp1 = GWSurrogates.NaturalCubicSpline(collect(ts), [u[1] for u in us])
        sp2 = GWSurrogates.NaturalCubicSpline(collect(ts), [u[2] for u in us])
        # Includes the last interval, where gwsurrogate's C++ spline evaluates the previous
        # interval's cubic for an isolated point instead.
        for x in (93.8, 94.5, 95.8, 96.5, 97.9, 98.0, 99.0, 100.0)
            v = GWSurrogates.interp4(ts..., us..., x)
            @test v ≈ SVector(GWSurrogates.spline_eval(sp1, x), GWSurrogates.spline_eval(sp2, x)) atol = 1.0e-13
        end
    end

    @testset "element_type promotes mixed tuples" begin
        D = ForwardDiff.Dual{Nothing, Float64, 1}
        d = D(1.0, ForwardDiff.Partials((1.0,)))
        # eltype of a tuple mixing duals and floats is an abstract union; element_type is concrete.
        @test !isconcretetype(eltype((d, 0.5, -0.3)))
        @test GWSurrogates.element_type((d, 0.5, -0.3)) === D
        @test GWSurrogates.element_type((1, 2.0)) === Float64
        @test GWSurrogates.element_type([1.0, 2.0]) === Float64
    end

    @testset "normalize_state" begin
        y = SVector(0.5, 1.0, -1.0, 2.0, 3.0, 0.3, 0.4, 0.0, -0.6, 0.0, 0.8)
        n = GWSurrogates.normalize_state(y, 0.25, 2.0)
        @test sum(abs2, n[1:4]) ≈ 1
        @test n[5] == 3.0
        @test sqrt(sum(abs2, n[6:8])) ≈ 0.25
        @test sqrt(sum(abs2, n[9:11])) ≈ 2.0
        z = GWSurrogates.normalize_state(SVector(1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 0.0, 0.0)
        @test z == SVector(1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    end
end

@testset "dynamics surrogate" begin
    ci = get(ENV, "CI", "false") == "true"
    for (model, params) in pairs(MODEL_PARAMS)
        fixture = joinpath(FIXTURE_DIR, "$(model)_dynamics.npz")
        if !isfile(fixture)
            msg = "dynamics reference fixture $fixture not found; generate it with python/generate_surrogates.py"
            ci ? error(msg) : @warn msg
            continue
        end
        @testset "$model" begin
            ref = npzread(fixture)
            sur = GWSurrogates.load_dynamics(model)
            @test sur.params == params
            @test sur.t == ref["t_ds"]
            @test length(GWSurrogates.output_times(sur)) == length(sur.t) - 3
            powers = zeros(GWSurrogates.n_fit_powers(sur.settings))

            @testset "right-hand side" begin
                q = ref["probe_q"]
                y = SVector{11}(ref["probe_y"])
                for (k, node) in enumerate(ref["probe_nodes"])
                    d = GWSurrogates.node_derivative(sur, node, q, y, powers)
                    @test maximum(abs, d - ref["probe_dydt_nodes"][k, :]) < 1.0e-14
                end
                for (k, t) in enumerate(ref["probe_times"])
                    d = GWSurrogates.time_derivative(sur, t, q, y, powers)
                    @test maximum(abs, d - ref["probe_dydt_times"][k, :]) < 1.0e-14
                end
            end

            @testset "reference time from omega" begin
                for (w, t_ref) in zip(ref["omega_refs"], ref["omega_t_refs"])
                    t = GWSurrogates.t_ref_from_omega(sur, w, ref["omega_q"], ref["omega_chiA"], ref["omega_chiB"])
                    @test t ≈ t_ref rtol = 1.0e-13
                end
                a, b = ref["omega_chiA"], ref["omega_chiB"]
                @test_throws ArgumentError GWSurrogates.t_ref_from_omega(sur, 0.3, 2.0, a, b)
                @test_throws ArgumentError GWSurrogates.t_ref_from_omega(sur, 0.001, 2.0, a, b)
            end

            @testset "trajectory $i" for i in 0:(ref["n_cases"] - 1)
                tr = ref["case$(i)_t_ref"]
                out = GWSurrogates.dynamics(
                    sur, ref["case$(i)_q"], ref["case$(i)_chiA0"], ref["case$(i)_chiB0"];
                    init_quat = ref["case$(i)_quat0"], init_orbphase = ref["case$(i)_orbphase0"],
                    t_ref = isnan(tr) ? nothing : tr,
                )
                @test maximum(abs, out.quat - ref["case$(i)_quat"]') < 1.0e-12
                @test maximum(abs, out.orbphase - ref["case$(i)_orbphase"]) < 1.0e-11
                @test maximum(abs, out.chiA - ref["case$(i)_chiA"]) < 1.0e-12
                @test maximum(abs, out.chiB - ref["case$(i)_chiB"]) < 1.0e-12
                # The per-step renormalization preserves the invariants.
                @test all(r -> sum(abs2, out.quat[r, :]) ≈ 1, axes(out.quat, 1))
                @test all(r -> sqrt(sum(abs2, out.chiA[r, :])) ≈ sqrt(sum(abs2, ref["case$(i)_chiA0"])), axes(out.chiA, 1))
            end

            # The kernels must be free of runtime dispatch for the AD backends.
            if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
                @testset "JET optimisation" begin
                    y = SVector{11}(ref["probe_y"])
                    a = (0.3, 0.2, 0.4)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.node_derivative(sur, 5, 2.0, y, powers)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.time_derivative(sur, -100.0, 2.0, y, powers)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.dynamics(sur, 2.0, a, a)
                    JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.dynamics(sur, 2.0, a, a; t_ref = -100.0)
                end
            end

            @testset "argument checks" begin
                a = (0.3, 0.2, 0.4)
                @test_throws ArgumentError GWSurrogates.dynamics(sur, 2.0, (0.0, 0.0, 1.5), a)
                @test_throws ArgumentError GWSurrogates.dynamics(sur, 2.0, a, a; t_ref = 500.0)
                @test_throws ArgumentError GWSurrogates.dynamics(sur, 2.0, a, a; t_ref = -5000.0)
            end

            @testset "ForwardDiff agrees with finite differences" begin
                fdm = central_fdm(5, 1)
                function objective(p; t_ref)
                    out = GWSurrogates.dynamics(
                        sur, p[1], (p[2], p[3], p[4]), (p[5], p[6], p[7]); init_orbphase = p[8], t_ref = t_ref,
                    )
                    return sum(out.quat[end, :] .* (1, 2, 3, 4)) + out.orbphase[end] + sum(out.chiA[end, :] .* (1, -1, 2)) +
                        sum(out.chiB[end, :] .* (2, 1, -1)) + sum(out.chiA[end ÷ 2, :])
                end
                p0 = [2.0, 0.3, 0.2, 0.4, -0.1, 0.5, -0.3, 0.4]
                for t_ref in (nothing, -100.0, -4000.0)
                    f(p) = objective(p; t_ref)
                    g = ForwardDiff.gradient(f, p0)
                    @test all(isfinite, g)
                    @test g ≈ grad(fdm, f, p0)[1] rtol = 1.0e-6 atol = 1.0e-7
                end
            end
        end
    end
end
