using Distributions
using QuadGK
using Random

@testset "UniformCosine" begin
    d = UniformCosine()
    total, _ = quadgk(x -> exp(logpdf(d, x)), 0.0, π)
    @test total ≈ 1.0 rtol = 1.0e-10
    for q in (0.0, 0.2, 0.5, 0.9, 1.0)
        @test cdf(d, quantile(d, q)) ≈ q atol = 1.0e-12
    end
    @test cdf(d, -1.0) == 0.0
    @test cdf(d, 4.0) ≈ 1.0
    @test logpdf(d, -0.1) == -Inf && logpdf(d, 3.2) == -Inf
    samples = rand(MersenneTwister(5), d, 100_000)
    @test all(s -> insupport(d, s), samples)
    @test mean(samples) ≈ π / 2 atol = 1.0e-2
end
