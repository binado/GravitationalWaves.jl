using Test, GWSEOBNR, ForwardDiff, Mooncake, ChainRulesCore, LinearAlgebra
const S = GWSEOBNR

function gradients(f, x)
    cache = Mooncake.prepare_gradient_cache(f, x)
    _, (_, g) = Mooncake.value_and_gradient!!(cache, f, x)
    return g, ForwardDiff.gradient(f, x)
end
relative(a, b) = norm(a - b) / norm(b)

@testset "Mooncake rules" begin
    @test Base.get_extension(S, :GWSEOBNRMooncakeExt) !== nothing

    @testset "real_fft" begin
        for n in (1, 2, 7, 8, 15, 16)
            x = [sin(0.7j) + 0.2cos(1.3j) for j in 1:n]
            w = [complex(cos(0.3k), sin(0.9k)) for k in 1:(div(n, 2) + 1)]
            loss = y -> sum(real(conj(w[k]) * z) + 0.13abs2(z) for (k, z) in enumerate(S.real_fft(y)))
            mooncake, forward = gradients(loss, x)
            @test mooncake ≈ forward rtol = 1.0e-12 atol = 1.0e-12
        end
    end

    @testset "implicit timing root" begin
        mooncake, forward = gradients(x -> S.continuous_root(t -> t^2 - x[1], 0.1, 3.0, t -> 2t), [2.0])
        @test mooncake ≈ forward rtol = 1.0e-12
    end

    @testset "solver-status validation keeps its errors" begin
        @test_throws S.SolveError S.check_nqc_solution([NaN], [1.0])
    end

    @testset "spectrum" begin
        grid = 0:0.125:1024
        generator = SEOBNRv5HM()
        function loss(x)
            h = S.evaluated_spectrum(
                grid, (;
                    mass_1 = x[1], mass_2 = x[2], chi_1 = x[3], chi_2 = x[4],
                    luminosity_distance = x[5], inclination = x[6], phase = x[7],
                ), generator,
            )
            return 1.0e22 * sum(real(h.plus[k]) + 0.3imag(h.plus[k]) + 0.7real(h.cross[k]) for k in 161:37:4001)
        end
        mooncake, forward = gradients(loss, [40.0, 20.0, 0.4, -0.3, 500.0, 0.7, 0.3])
        @test relative(mooncake, forward) < 1.0e-8
    end
end
