using Reactant
using Enzyme: Enzyme

Reactant.set_default_backend("cpu")

_host(x::Number) = Reactant.to_number(x)
_host(x::Union{Tuple, NamedTuple}) = map(_host, x)
_host(x) = x

function _compiled_close(actual, expected; rtol = 1.0e-12, atol = 1.0e-13)
    return if expected isa Union{Tuple, NamedTuple}
        @test keys(actual) == keys(expected)
        foreach((a, b) -> _compiled_close(a, b; rtol, atol), values(actual), values(expected))
    elseif expected isa Bool || !(expected isa Number)
        @test actual == expected
    else
        @test isapprox(actual, expected; rtol, atol)
    end
end

# Each executable must handle multiple tracked values of the same structure.
function _test_compiled(f, samples; kwargs...)
    tracked = Reactant.to_rarray(first(samples); track_numbers = true)
    executable = Reactant.@compile f(tracked...)
    for sample in samples
        args = Reactant.to_rarray(sample; track_numbers = true)
        _compiled_close(_host(executable(args...)), f(sample...); kwargs...)
    end
    return executable
end

function _test_compiled_gradient(f, samples; rtol = 1.0e-10, atol = 1.0e-11)
    # A vector permits comparison of every partial against an independent AD backend.
    tracked = Reactant.to_rarray(collect(first(samples)))
    objective(v) = Reactant.@allowscalar sum(_flatten_numbers(f(v...)))
    executable = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, tracked)
    for sample in samples
        v = collect(sample)
        actual = Array(only(executable(Enzyme.Reverse, objective, Reactant.to_rarray(v))))
        @test isapprox(actual, ForwardDiff.gradient(objective, v); rtol, atol)
    end
    return
end

_flatten_numbers(x::Number) = (x,)
_flatten_numbers(x::Union{Tuple, NamedTuple}) = mapreduce(_flatten_numbers, (a, b) -> (a..., b...), values(x))

@testset "Reactant CPU" begin
    @testset "public numerical functions and gradients" begin
        m1, m2 = 36.0, 29.0
        n1, n2 = 42.0, 23.0
        mc, nc = chirp_mass(m1, m2), chirp_mass(n1, n2)
        cases = (
            (
                (
                    chirp_mass, total_mass, mass_ratio, symmetric_mass_ratio, reduced_mass,
                    component_masses_to_chirp_mass_and_mass_ratio,
                    component_masses_to_total_mass_and_mass_ratio,
                ), ((m1, m2), (n1, n2)),
            ),
            ((total_mass_and_mass_ratio_to_component_masses,), ((65.0, m2 / m1), (70.0, 0.6))),
            (
                (chirp_mass_and_mass_ratio_to_component_masses, chirp_mass_and_mass_ratio_to_total_mass),
                ((mc, m2 / m1), (nc, n2 / n1)),
            ),
            ((chirp_mass_and_total_mass_to_symmetric_mass_ratio,), ((mc, 65.0), (nc, 65.0))),
            ((symmetric_mass_ratio_to_mass_ratio,), ((0.2,), (0.24,))),
            ((mass_1_and_chirp_mass_to_mass_ratio,), ((m1, mc), (n1, nc))),
            ((mass_2_and_chirp_mass_to_mass_ratio,), ((m2, mc), (n2, nc))),
            ((source_to_detector_mass, detector_to_source_mass), ((30.0, 0.5), (42.0, 0.7))),
            ((effective_spin,), ((0.3, -0.2, 0.8), (-0.4, 0.6, 0.5))),
            ((effective_precessing_spin,), ((0.3, 0.4, 0.8), (0.1, 0.6, 0.5))),
            ((effective_precessing_spin,), ((0.3, 0.4, 0.2, 0.1, 0.8), (0.1, 0.2, 0.4, 0.5, 0.7))),
            ((aligned_spins_to_lal_frame,), ((0.4, 0.3, 0.0, 0.2, Float64(π)), (0.8, 0.5, Float64(π), 0.4, 0.0))),
            (
                (transform_precessing_initial_conditions,),
                (
                    (0.4, 2.5, 0.7, 1.9, 1.1, 0.6, 0.4, 36.0, 29.0, 50.0, 1.3),
                    (0.8, 1.5, 0.9, 1.4, 0.7, 0.5, 0.3, 42.0, 23.0, 40.0, 0.9),
                ),
            ),
        )
        for (functions, samples) in cases, f in functions
            @testset "$f/$(length(first(samples)))" begin
                _test_compiled(f, samples)
                _test_compiled_gradient(f, samples)
            end
        end
        _test_compiled(is_aligned, ((0.5, 0.7), (0.0, 0.7), (0.5, 0.0), (0.5, Float64(π))))
        # Boundary behavior is tested without derivatives at singularities/ties.
        _test_compiled(symmetric_mass_ratio_to_mass_ratio, ((0.2,), (0.25,)))
        _test_compiled(effective_precessing_spin, ((0.3, 0.4, 0.8), (0.0, 0.0, 1.0)))
        _test_compiled(effective_precessing_spin, ((0.3, 0.4, 0.2, 0.1, 0.8), (0.0, 0.0, 0.0, 0.0, 1.0)))
    end

    @testset "conversion registry" begin
        function parameters(m1, m2, a1, a2, t1, t2, z)
            p = (; mass_1 = m1, mass_2 = m2, a_1 = a1, a_2 = a2, tilt_1 = t1, tilt_2 = t2, redshift = z)
            return convert_parameters(p, (GWParameterConversion.MASS_KEYS..., :chi_eff, :chi_p))
        end
        full = (
            parameters(36.0, 29.0, 0.6, 0.4, 0.7, 1.9, 0.5),
            parameters(42.0, 23.0, 0.5, 0.3, 0.9, 1.4, 0.7),
        )
        for c in GWParameterConversion.REGISTRY
            f(nt) = convert_parameters(nt, c.gives; check_consistency = false)
            samples = map(p -> (NamedTuple{c.needs}(map(k -> getproperty(p, k), c.needs)),), full)
            _test_compiled(f, samples)
            # Check the planner against the registry entry, not just host planning.
            for (nt,) in samples
                expected = GWParameterConversion._tuple(c.f(values(nt)...))
                output = f(nt)
                _compiled_close(map(k -> getproperty(output, k), c.gives), expected)
            end
            g(xs...) = f(NamedTuple{c.needs}(xs))
            _test_compiled_gradient(g, map(x -> values(only(x)), samples))
        end
        targets = (:mass_1, :mass_2, :chirp_mass_source, :total_mass_source, :chi_eff, :chi_p)
        chain(nt) = convert_parameters(nt, targets; check_consistency = false)
        inputs = map(full) do p
            (;
                chirp_mass_source = p.chirp_mass_source, mass_ratio = p.mass_ratio,
                redshift = p.redshift, a_1 = p.a_1, a_2 = p.a_2, tilt_1 = p.tilt_1, tilt_2 = p.tilt_2,
                unrelated = "kept",
            )
        end
        _test_compiled(chain, map(x -> (x,), inputs))
        chain_gradient(mc, q, z, a1, a2, t1, t2) = chain(
            (;
                chirp_mass_source = mc,
                mass_ratio = q, redshift = z, a_1 = a1, a_2 = a2, tilt_1 = t1, tilt_2 = t2,
            )
        )
        _test_compiled_gradient(
            chain_gradient, map(
                p -> (
                    p.chirp_mass_source, p.mass_ratio,
                    p.redshift, p.a_1, p.a_2, p.tilt_1, p.tilt_2,
                ), full
            )
        )
        # Consistent redundant fields must remain supplied, including their order.
        redundant(nt) = convert_parameters(nt, (:total_mass, :chi_eff); check_consistency = false)
        samples = map(
            p -> (
                (;
                    mass_1 = p.mass_1, mass_2 = p.mass_2, mass_ratio = p.mass_ratio,
                    a_1 = p.a_1, a_2 = p.a_2, tilt_1 = p.tilt_1, tilt_2 = p.tilt_2,
                ),
            ), full
        )
        _test_compiled(redundant, samples)
        for (nt,) in samples
            @test redundant(nt).mass_ratio === nt.mass_ratio
        end
    end

    @testset "spin frames across branches" begin
        base = (; mass_1 = 36.0, mass_2 = 29.0, theta_jn = 0.4, phase = 1.3)
        p = merge(base, (; a_1 = 0.6, a_2 = 0.4, tilt_1 = 0.7, tilt_2 = 1.9, phi_12 = 1.1, phi_jl = 2.5))
        inputs = (
            p, merge(p, (; tilt_1 = 0.0, tilt_2 = 0.0)),
            merge(p, (; tilt_1 = Float64(π), tilt_2 = Float64(π))),
            merge(p, (; a_1 = 0.0, a_2 = 0.0)),
        )
        executable = _test_compiled(spins_to_lal_frame, map(x -> (x, 50.0), inputs))
        for x in inputs[2:end]
            s = _host(
                executable(
                    Reactant.to_rarray(x; track_numbers = true),
                    Reactant.to_rarray(50.0; track_numbers = true)
                )
            )
            @test (s.s1x, s.s1y, s.s2x, s.s2y) == (0.0, 0.0, 0.0, 0.0)
            @test s.iota == x.theta_jn
        end
        # The defaults are static keys; their numerical results remain dynamic.
        _test_compiled(spins_to_lal_frame, ((base, 50.0), (merge(base, (; theta_jn = 0.8)), 40.0)))
        names = keys(p)
        frame(xs...) = spins_to_lal_frame(NamedTuple{names}(xs), 50.0)
        _test_compiled_gradient(frame, (values(p), values(merge(p, (; tilt_1 = 0.9, theta_jn = 0.8)))))

        # Opposite equal in-plane spins cancel in J, so J = L = N.
        # This stays inside the fallback threshold without inferring an angle
        # from rounded host values near the threshold.
        fallback_p = merge(
            p, (;
                mass_1 = 30.0, mass_2 = 30.0, a_1 = 0.5, a_2 = 0.5,
                tilt_1 = Float64(π) / 2, tilt_2 = Float64(π) / 2, phi_12 = Float64(π), theta_jn = 0.0,
            )
        )
        fallback = (
            fallback_p.theta_jn, fallback_p.phi_jl, fallback_p.tilt_1,
            fallback_p.tilt_2, fallback_p.phi_12, fallback_p.a_1, fallback_p.a_2,
            fallback_p.mass_1, fallback_p.mass_2, 50.0, fallback_p.phase,
        )
        _compiled_close(
            _host(
                executable(
                    Reactant.to_rarray(fallback_p; track_numbers = true),
                    Reactant.to_rarray(50.0; track_numbers = true)
                )
            ), spins_to_lal_frame(fallback_p, 50.0)
        )
        _test_compiled(
            transform_precessing_initial_conditions,
            (fallback, (0.3, fallback[2:end]...))
        )
        fallback_result = transform_precessing_initial_conditions(fallback...)
        @test fallback_result.iota == 0.0
        # The fixed y × L direction sets the azimuth when the sightline cannot.
        @test hypot(fallback_result.s1x, fallback_result.s1y) ≈ 0.5
        azimuth = p.phi_jl + Float64(π) / 2 - p.phase
        @test fallback_result.s1x ≈ 0.5 * cos(azimuth) rtol = 1.0e-12
        @test fallback_result.s1y ≈ 0.5 * sin(azimuth) rtol = 1.0e-12

        direct = (0.0, 1.0, 0.0, 0.0, 1.1, 0.6, 0.4, 36.0, 29.0, 50.0, 1.3)
        _test_compiled(transform_precessing_initial_conditions, (direct, (Float64(π), direct[2:end]...)))
        host = transform_precessing_initial_conditions(direct...)
        @test all(isfinite, values(host))
        @test host.s1z ≈ 0.6
        @test host.s2z ≈ 0.4
    end

    @testset "coordinates, Jacobians and inverses" begin
        for (forward, key) in (
                (chirp_mass_and_mass_ratio_to_component_masses, :chirp_mass),
                (total_mass_and_mass_ratio_to_component_masses, :total_mass),
            )
            backward = inverse(forward)
            xs = (
                NamedTuple{(key, :mass_ratio)}((28.0, 0.8)),
                NamedTuple{(key, :mass_ratio)}((35.0, 0.6)),
                NamedTuple{(key, :mass_ratio)}((30.0, 1.0)),
            )
            for f in (forward, backward)
                inputs = f === forward ? xs : map(forward, xs)
                for samples in (inputs, map(x -> NamedTuple{reverse(keys(x))}(reverse(values(x))), inputs))
                    _test_compiled(f, map(x -> (x,), samples))
                    jac(x) = with_logabsdet_jacobian(f, x)
                    _test_compiled(jac, map(x -> (x,), samples))
                    roundtrip(x) = inverse(f)(f(x))
                    _test_compiled(roundtrip, map(x -> (x,), samples))
                    function jacobian_roundtrip(x)
                        y, j = with_logabsdet_jacobian(f, x)
                        restored, invj = with_logabsdet_jacobian(inverse(f), y)
                        return restored, j + invj
                    end
                    _test_compiled(jacobian_roundtrip, map(x -> (x,), samples))
                    names = keys(first(samples))
                    objective(v...) = with_logabsdet_jacobian(f, NamedTuple{names}(v))
                    _test_compiled_gradient(objective, map(values, samples[1:2]))
                end
            end
        end
        f = chirp_mass_and_mass_ratio_to_component_masses
        g = component_masses_to_total_mass_and_mass_ratio
        composed(x) = with_logabsdet_jacobian(g ∘ f, x)
        _test_compiled(
            composed, (
                ((; chirp_mass = 28.0, mass_ratio = 0.8),),
                ((; chirp_mass = 35.0, mass_ratio = 0.6),),
            )
        )
    end
end
