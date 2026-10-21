using FiniteDifferences
using ForwardDiff

@testset "guarded norms" begin
    fdm = central_fdm(5, 1)

    @testset "values" begin
        @test SurrogateWaveforms.guarded_sqrt(4.0) == 2.0
        @test SurrogateWaveforms.guarded_sqrt(0.0) == 0.0
        @test SurrogateWaveforms.guarded_sqrt(-1.0) == 0.0
        v = SurrogateWaveforms.rescale_to_norm((3.0, 0.0, 4.0), 0.5)
        @test all(v .≈ (0.3, 0.0, 0.4))
        @test SurrogateWaveforms.rescale_to_norm((0.0, 0.0, 0.0), 0.5) == (0.0, 0.0, 0.0)
    end

    @testset "zero spin gives finite gradients" begin
        norm_objective(c) = SurrogateWaveforms.guarded_sqrt(sum(abs2, c))
        scale_objective(c) = sum(SurrogateWaveforms.rescale_to_norm((c[1], c[2], c[3]), 0.7))
        # The norm has zero slope at the guard; the rescale is the identity on a zero vector.
        @test ForwardDiff.gradient(norm_objective, zeros(3)) == zeros(3)
        @test ForwardDiff.gradient(scale_objective, zeros(3)) == ones(3)
    end

    @testset "matches finite differences away from zero" begin
        c = [0.3, -0.4, 0.2]
        norm_objective(c) = SurrogateWaveforms.guarded_sqrt(sum(abs2, c))
        scale_objective(c) = sum(abs2, SurrogateWaveforms.rescale_to_norm((c[1], c[2], c[3]), 0.7) .* (1.0, 2.0, 3.0))
        for f in (norm_objective, scale_objective)
            @test ForwardDiff.gradient(f, c) ≈ grad(fdm, f, c)[1] atol = 1.0e-9
        end
    end
end
