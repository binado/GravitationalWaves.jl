using Distributions
using Random

@testset "InterpolatedDistribution table" begin
    x = [0.0, 1.0, 2.0, 4.0]
    y = [0.0, 2.0, 2.0, 0.0]
    d = InterpolatedDistribution(x, y)
    @test normalizer(d) == GWDistributions._trapz(y, x) == 5.0
    @test d.cdf_grid[1] == 0.0
    @test last(d.cdf_grid) === 1.0
    @test d.cdf_grid ≈ [0.0, 1 / 5, 3 / 5, 1.0]
    @test quantile(d, 1.0) == x[end]
    @test quantile(d, 0.0) == x[1]
    @test minimum(d) == 0.0 && maximum(d) == 4.0
    @test insupport(d, 4.0) && !insupport(d, 4.1)

    @test logpdf(d, 1.5) ≈ log(2 / 5)
    @test logpdf(d, 0.5) ≈ log(1 / 5)
    # Off the grid, and where the table is zero, the log density is -Inf.
    @test logpdf(d, -0.1) == -Inf
    @test logpdf(d, 4.1) == -Inf
    @test logpdf(d, 0.0) == -Inf
    @test logpdf(d, 4.0) == -Inf
end

@testset "InterpolatedDistribution cdf / quantile / rand" begin
    x = collect(range(0.0, 3.0; length = 31))
    y = @. 1 + x^2
    d = InterpolatedDistribution(x, y)
    for q in (0.0, 0.05, 0.3, 0.5, 0.97, 1.0)
        @test cdf(d, quantile(d, q)) ≈ q atol = 1.0e-12
    end
    @test cdf(d, -1.0) == 0.0
    @test cdf(d, 4.0) == 1.0
    samples = rand(MersenneTwister(7), d, 200)
    @test all(s -> insupport(d, s), samples)
    # The linear interpolant of y is integrated exactly by the trapezoid CDF.
    @test cdf(d, 1.5) ≈ GWDistributions._trapz(y[1:16], x[1:16]) / normalizer(d)

    @test_throws DimensionMismatch InterpolatedDistribution([0.0, 1.0], [1.0])
    @test_throws ArgumentError InterpolatedDistribution([0.0], [1.0])
end
