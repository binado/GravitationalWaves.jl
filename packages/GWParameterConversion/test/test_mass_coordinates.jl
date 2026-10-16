using ChangesOfVariables: with_logabsdet_jacobian
using InverseFunctions: inverse
using LinearAlgebra: logabsdet

@testset "mass coordinates" begin
    transforms = (
        (chirp_mass_and_mass_ratio_to_component_masses, component_masses_to_chirp_mass_and_mass_ratio, :chirp_mass),
        (total_mass_and_mass_ratio_to_component_masses, component_masses_to_total_mass_and_mass_ratio, :total_mass),
    )
    components = [(36.0, 29.0), (30.0, 30.0), (100.0, 0.1), (0.01, 0.002)]
    close(a, b) = keys(a) == keys(b) && all(isapprox.(values(a), values(b); rtol = 1.0e-12))
    # NamedTuple field order defines the numerical input axes for this oracle.
    jacobian(f, p) = ForwardDiff.jacobian(
        x -> collect(values(f(NamedTuple{keys(p)}(Tuple(x))))), collect(values(p))
    )

    @testset "$mass_key" for (forward, backward, mass_key) in transforms
        @test inverse(forward) === backward
        @test inverse(backward) === forward
        @test inverse(inverse(forward)) === forward

        for (m1, m2) in components
            p = (; mass_1 = m1, mass_2 = m2)
            x = @inferred backward(p)
            @test keys(x) == (mass_key, :mass_ratio)
            @test values(x) == backward(m1, m2)
            @test values(@inferred forward(x)) == forward(values(x)...)
            @test close(forward(x), p)
            @test close(backward(forward(x)), x)

            for (f, input) in ((forward, x), (backward, p))
                y, ladj = @inferred with_logabsdet_jacobian(f, input)
                @test y == f(input)
                @test ladj ≈ logabsdet(jacobian(f, input))[1] rtol = 1.0e-12
                restored, inverse_ladj = with_logabsdet_jacobian(inverse(f), y)
                @test close(restored, input)
                @test ladj ≈ -inverse_ladj rtol = 1.0e-12

                reversed = NamedTuple{reverse(keys(input))}(reverse(values(input)))
                @test f(reversed) == y
                @test with_logabsdet_jacobian(f, reversed) == (y, ladj)
                @test ladj ≈ logabsdet(jacobian(f, reversed))[1] rtol = 1.0e-12

                # Differentiation must work through both outputs of the combined call.
                z = collect(values(input))
                gradient = ForwardDiff.gradient(z) do v
                    nt = NamedTuple{keys(input)}(Tuple(v))
                    output, logjac = with_logabsdet_jacobian(f, nt)
                    return sum(values(output)) + logjac
                end
                h = 1.0e-5
                component_input = keys(input) == (:mass_1, :mass_2)
                objective(t) = begin
                    nt = NamedTuple{keys(input)}((z[1] * t, component_input ? z[2] * t : z[2]))
                    output, logjac = with_logabsdet_jacobian(f, nt)
                    sum(values(output)) + logjac
                end
                # Scale masses together so equal-mass inputs stay in the domain.
                fd = (objective(1 + h) - objective(1 - h)) / (2h)
                directional = gradient[1] * z[1] + (component_input ? gradient[2] * z[2] : 0)
                @test directional ≈ fd rtol = 1.0e-7
                @test all(isfinite, gradient)
            end
        end

        @testset "domain and key validation" begin
            valid = NamedTuple{(mass_key, :mass_ratio)}((28.0, 0.8))
            for bad in ((;), (; mass_ratio = 0.8), merge(valid, (; extra = 1)), (; wrong_mass = 28.0, mass_ratio = 0.8))
                @test_throws ArgumentError forward(bad)
                @test_throws ArgumentError with_logabsdet_jacobian(forward, bad)
            end
            @test_throws ArgumentError forward(NamedTuple{keys(valid)}(("28", 0.8)))
            @test_throws ArgumentError backward((; mass_1 = 36.0, mass_2 = 29.0, extra = 1))
            @test_throws ArgumentError backward((; mass_1 = 36.0))
            @test_throws ArgumentError backward((; mass_1 = 36.0, mass_2 = 29.0im))
            for m in (0.0, -1.0, Inf, -Inf, NaN)
                bad = NamedTuple{keys(valid)}((m, 0.8))
                @test_throws DomainError forward(bad)
                @test_throws DomainError with_logabsdet_jacobian(forward, bad)
                @test_throws DomainError backward((; mass_1 = m, mass_2 = 1.0))
                @test_throws DomainError backward((; mass_1 = 36.0, mass_2 = m))
            end
            for q in (0.0, -0.1, 1.1, Inf, -Inf, NaN)
                bad = NamedTuple{keys(valid)}((28.0, q))
                @test_throws DomainError forward(bad)
                @test_throws DomainError with_logabsdet_jacobian(forward, bad)
            end
            @test_throws DomainError backward((; mass_1 = 29.0, mass_2 = 36.0))
            @test_throws DomainError with_logabsdet_jacobian(backward, (; mass_1 = 29.0, mass_2 = 36.0))
        end

        # Retain precision and AD types, including in validation and logarithms.
        big = NamedTuple{(mass_key, :mass_ratio)}((big"28.0", big"0.8"))
        y, ladj = @inferred with_logabsdet_jacobian(forward, big)
        @test all(v -> v isa BigFloat, values(y))
        @test ladj isa BigFloat
        dual = NamedTuple{keys(big)}((ForwardDiff.Dual(28.0, 1.0), ForwardDiff.Dual(0.8, 0.0)))
        ydual, ldual = @inferred with_logabsdet_jacobian(forward, dual)
        @test all(v -> v isa ForwardDiff.Dual, values(ydual))
        @test ldual isa ForwardDiff.Dual

        @testset "AD at the equal-mass boundary" begin
            # Nested differentiation must validate primal coordinates, including
            # negative seeds and unequal tangents at equal component masses.
            for q in (ForwardDiff.Dual(1.0, -1.0), ForwardDiff.Dual(1.0, 1.0))
                x = NamedTuple{(mass_key, :mass_ratio)}((28.0, q))
                @test all(isfinite, values(forward(x)))
            end
            p = (; mass_1 = ForwardDiff.Dual(30.0, -1.0), mass_2 = ForwardDiff.Dual(30.0, 1.0))
            @test all(isfinite, values(backward(p)))
            hessian = ForwardDiff.hessian([28.0, 1.0]) do v
                x = NamedTuple{(mass_key, :mass_ratio)}(Tuple(v))
                y, ladj = with_logabsdet_jacobian(forward, x)
                return sum(values(y)) + ladj
            end
            @test all(isfinite, hessian)
            @test_throws DomainError forward(NamedTuple{(mass_key, :mass_ratio)}((28.0, ForwardDiff.Dual(0.0, 1.0))))
        end
    end

    @testset "composition" begin
        f = chirp_mass_and_mass_ratio_to_component_masses
        g = component_masses_to_total_mass_and_mass_ratio
        x = (; chirp_mass = 28.0, mass_ratio = 0.8)
        y, j1 = with_logabsdet_jacobian(f, x)
        z, j2 = with_logabsdet_jacobian(g, y)
        combined, ladj = @inferred with_logabsdet_jacobian(g ∘ f, x)
        @test combined == z
        @test z.total_mass ≈ chirp_mass_and_mass_ratio_to_total_mass(x.chirp_mass, x.mass_ratio)
        @test ladj ≈ j1 + j2
        @test ladj ≈ logabsdet(jacobian(g ∘ f, x))[1]
        @test close(inverse(g ∘ f)(combined), x)
    end

    @testset "independent extension loading" begin
        # Fresh processes exclude optional imports from this test session.
        for deps in ((), (:ChangesOfVariables,), (:InverseFunctions,), (:ForwardDiff,), (:ChangesOfVariables, :InverseFunctions), (:InverseFunctions, :ChangesOfVariables), (:ForwardDiff, :ChangesOfVariables, :InverseFunctions))
            imports = join(("using $dep" for dep in deps), "\n")
            cov = :ChangesOfVariables in deps
            inv = :InverseFunctions in deps
            fd = :ForwardDiff in deps
            script = """
            using GWParameterConversion
            @assert Base.get_extension(GWParameterConversion, :GWParameterConversionChangesOfVariablesExt) === nothing
            @assert Base.get_extension(GWParameterConversion, :GWParameterConversionInverseFunctionsExt) === nothing
            @assert Base.get_extension(GWParameterConversion, :GWParameterConversionForwardDiffExt) === nothing
            $imports
            @assert (Base.get_extension(GWParameterConversion, :GWParameterConversionChangesOfVariablesExt) !== nothing) == $cov
            @assert (Base.get_extension(GWParameterConversion, :GWParameterConversionInverseFunctionsExt) !== nothing) == $inv
            @assert (Base.get_extension(GWParameterConversion, :GWParameterConversionForwardDiffExt) !== nothing) == $fd
            f = chirp_mass_and_mass_ratio_to_component_masses
            p = (; chirp_mass = 28.0, mass_ratio = 0.8)
            @assert keys(f(p)) == (:mass_1, :mass_2)
            """
            if cov
                script *= "\n@assert ChangesOfVariables.with_logabsdet_jacobian(f, p)[1] == f(p)"
            end
            if inv
                script *= "\n@assert InverseFunctions.inverse(f) === component_masses_to_chirp_mass_and_mass_ratio"
            end
            if fd
                script *= "\n@assert isfinite(ForwardDiff.derivative(q -> f((; chirp_mass = 28.0, mass_ratio = q)).mass_1, 1.0))"
            end
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`
            @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
        end
    end
end
