using TestItems

@testsnippet Setup_delay begin
    using Distributions
    using ForwardDiff
    using QuadGK
    using Random
end

@testitem "PowerLawDelay normalization and round trip" setup = [Common, Setup_delay] begin
    for α in (-2.0, -1.0, -0.3, 1.5)
        d = PowerLawDelay(α, 0.02, 13.0)
        total, _ = quadgk(t -> exp(logpdf(d, t)), 0.02, 13.0; rtol = 1.0e-10)
        @test total ≈ 1.0 rtol = 1.0e-8
        for q in (0.0, 0.1, 0.5, 0.9, 1.0)
            @test cdf(d, quantile(d, q)) ≈ q atol = 1.0e-12
        end
        @test quantile(d, 0.0) ≈ 0.02
        @test quantile(d, 1.0) ≈ 13.0
        @test cdf(d, 0.001) == 0.0
        @test cdf(d, 20.0) == 1.0
        @test logpdf(d, 0.01) == -Inf
        @test logpdf(d, 14.0) == -Inf
        @test all(s -> insupport(d, s), rand(MersenneTwister(1), d, 50))
    end
    @test_throws ArgumentError PowerLawDelay(-1.0, 0.0, 13.0)
    @test_throws ArgumentError PowerLawDelay(-1.0, 13.0, 0.02)
end

@testitem "PowerLawDelay matches the closed forms away from α = -1" setup = [Common, Setup_delay] begin
    low, high = 0.02, 13.0
    for α in (-2.5, -0.5, 0.7)
        d = PowerLawDelay(α, low, high)
        β = 1 + α
        for t in (0.05, 1.0, 10.0)
            @test cdf(d, t) ≈ (t^β - low^β) / (high^β - low^β) rtol = 1.0e-12
            @test logpdf(d, t) ≈ log(β / (high^β - low^β)) + α * log(t) rtol = 1.0e-12
        end
        for q in (0.1, 0.5, 0.9)
            @test quantile(d, q) ≈ (low^β + q * (high^β - low^β))^(1 / β) rtol = 1.0e-12
        end
    end
    d = PowerLawDelay(-1.0, low, high)
    @test cdf(d, 1.0) ≈ log(1.0 / low) / log(high / low)
    @test quantile(d, 0.5) ≈ low * (high / low)^0.5
end

@testitem "PowerLawDelay is smooth through α = -1 under AD" setup = [Common, Setup_delay] begin
    t, q = 1.0, 0.4
    dcdf = α -> ForwardDiff.derivative(a -> cdf(PowerLawDelay(a, 0.02, 13.0), t), α)
    dq = α -> ForwardDiff.derivative(a -> quantile(PowerLawDelay(a, 0.02, 13.0), q), α)
    for f in (dcdf, dq)
        centre = f(-1.0)
        @test isfinite(centre)
        @test f(-1.0 - 1.0e-4) ≈ centre rtol = 1.0e-3
        @test f(-1.0 + 1.0e-4) ≈ centre rtol = 1.0e-3
        # Right next to the removable singularity: still finite and close.
        @test f(-1.0 + 1.0e-9) ≈ centre rtol = 1.0e-3
    end
    # And the derivative agrees with a central finite difference.
    fd = (cdf(PowerLawDelay(-1.0 + 1.0e-6, 0.02, 13.0), t) - cdf(PowerLawDelay(-1.0 - 1.0e-6, 0.02, 13.0), t)) / 2.0e-6
    @test dcdf(-1.0) ≈ fd rtol = 1.0e-6
end
