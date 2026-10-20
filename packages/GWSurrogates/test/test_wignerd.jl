import Quaternionic
import SphericalFunctions
using FiniteDifferences
using ForwardDiff
using LinearAlgebra
using Random

@testset "Wigner-D" begin
    oracle(R, l) = SphericalFunctions.D(R, l)[l]

    function rotor_cases()
        rng = MersenneTwister(7)
        cases = [
            ("identity", (1.0, 0.0, 0.0, 0.0)),
            ("z rotation", (cos(0.4), 0.0, 0.0, sin(0.4))),
            ("beta = pi", (0.0, cos(0.3), sin(0.3), 0.0)),
            ("x rotation", (cos(0.9), sin(0.9), 0.0, 0.0)),
        ]
        for k in 1:5
            q = randn(rng, 4)
            push!(cases, ("random $k", Tuple(q ./ norm(q))))
        end
        return cases
    end

    @testset "matches SphericalFunctions.D, l = $l" for l in 0:5
        for (name, (w, x, y, z)) in rotor_cases()
            R = Quaternionic.Quaternion(w, x, y, z)
            ref = [oracle(R, l)[mp, m] for mp in (-l):l, m in (-l):l]
            @test maximum(abs, GWSurrogates.wigner_D(l, w, x, y, z) - ref) < 1.0e-13
        end
    end

    @testset "is unitary" begin
        w, x, y, z = (0.5, -0.5, 0.5, 0.5)
        D = GWSurrogates.wigner_D(4, w, x, y, z)
        @test D' * D ≈ I atol = 1.0e-13
    end

    @testset "argument checks" begin
        @test_throws ArgumentError GWSurrogates.wigner_D(6, 1.0, 0.0, 0.0, 0.0)
        @test_throws ArgumentError GWSurrogates.wigner_D_element(2, 3, 0, 1.0 + 0im, 0.0im)
    end

    @testset "ForwardDiff agrees with finite differences" begin
        fdm = central_fdm(5, 1)
        # Real scalar objective over a fixed weighted sum of all blocks, at a normalised quaternion.
        rng = MersenneTwister(11)
        weights = [randn(rng, ComplexF64, 2l + 1, 2l + 1) for l in 2:5]
        function objective(q)
            n = sqrt(sum(abs2, q))
            w, x, y, z = (q[1] / n, q[2] / n, q[3] / n, q[4] / n)
            sum(real(sum(weights[l - 1] .* GWSurrogates.wigner_D(l, w, x, y, z))) for l in 2:5)
        end
        points = Dict(
            "generic" => [0.3, -0.2, 0.7, 0.4],
            "identity" => [1.0, 0.0, 0.0, 0.0],
            "z rotation" => [cos(0.4), 0.0, 0.0, sin(0.4)],
            "beta = pi" => [0.0, cos(0.3), sin(0.3), 0.0],
        )
        for (name, q) in points
            @testset "$name" begin
                g = ForwardDiff.gradient(objective, q)
                @test all(isfinite, g)
                @test g ≈ grad(fdm, objective, q)[1] atol = 1.0e-8
            end
        end
    end
end
