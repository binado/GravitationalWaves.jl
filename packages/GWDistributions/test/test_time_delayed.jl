using BackgroundCosmology
using Distributions
using ForwardDiff
using QuadGK
using Random

const _TD_PARAMS = (H0 = 67.66, Ωm = 0.30966, γ = 2.7, κ = 3.0, z_peak = 2.0, local_merger_rate = 17.0)
const _TD_DELAY = PowerLawDelay(-1.0, 0.02, 13.0)
_td_ψ(z, p) = GWDistributions._madau_dickinson_source_frame(z, p)

@testset "_redshift_at_lookback_time inverts lookback_time" begin
    for c in (
            cosmology((H0 = 67.66, Ωm = 0.30966)),
            cosmology((H0 = 70.0, Ωm = 0.3, w0 = -0.9)),
            cosmology((H0 = 70.0, Ωm = 0.3, w0 = -0.9, wa = 0.2)),
        )
        for z in (0.0, 0.1, 1.0, 5.0, 19.0)
            t = lookback_time(c, z)
            @test GWDistributions._redshift_at_lookback_time(c, t) ≈ z atol = 1.0e-8 rtol = 1.0e-8
        end
        @test GWDistributions._redshift_at_lookback_time(c, age(c, 0) + 1) == Inf
    end
    # The closed form is differentiable in the cosmology.
    g = ForwardDiff.derivative(Ωm -> GWDistributions._redshift_at_lookback_time(cosmology((H0 = 70.0, Ωm = Ωm)), 5.0), 0.3)
    @test isfinite(g)
    # So is the Newton path.
    gw = ForwardDiff.derivative(Ωm -> GWDistributions._redshift_at_lookback_time(cosmology((H0 = 70.0, Ωm = Ωm, w0 = -0.9)), 5.0), 0.3)
    @test isfinite(gw)
    @test gw ≈ g rtol = 0.2
end

@testset "time_delayed_merger_rate" begin
    kw = (; formation_rate = _td_ψ, delay = _TD_DELAY)
    @test time_delayed_merger_rate(0.0, _TD_PARAMS; kw...) ≈ _td_ψ(0.0, _TD_PARAMS)

    # Reference: adaptive quadrature of the convolution in τ.
    c = cosmology(_TD_PARAMS)
    t_cut = lookback_time(c, 20.0)
    function convolution(z)
        t = lookback_time(c, z)
        upper = min(13.0, t_cut - t)
        val, _ = quadgk(0.02, upper; rtol = 1.0e-12) do τ
            exp(logpdf(_TD_DELAY, τ)) *
                _td_ψ(GWDistributions._redshift_at_lookback_time(c, t + τ), _TD_PARAMS)
        end
        return val
    end
    for z in (0.3, 1.0, 3.0)
        expected = convolution(z) / convolution(0.0) * _td_ψ(0.0, _TD_PARAMS)
        @test time_delayed_merger_rate(z, _TD_PARAMS; kw...) ≈ expected rtol = 1.0e-8
    end
end

@testset "TimeDelayedRedshiftDistribution" begin
    d = madau_dickinson_time_delayed_redshift_distribution(;
        params = _TD_PARAMS, time_delay_distribution = _TD_DELAY, maximum_redshift = 3.0, n_grid = 61
    )
    @test d isa TimeDelayedRedshiftDistribution
    @test d isa RedshiftDistribution
    @test merger_rate(d, 0.0) ≈ _TD_PARAMS.local_merger_rate
    @test merger_rate(d, 1.0) ≈ time_delayed_merger_rate(
        1.0, _TD_PARAMS; formation_rate = _td_ψ, delay = _TD_DELAY
    )
    # The accessor returns the formation rate, not the delayed rate.
    @test source_frame_distribution(d, 1.0, _TD_PARAMS) == _td_ψ(1.0, _TD_PARAMS)
    @test total_merger_rate(d) > 0
    @test cdf(d, quantile(d, 0.6)) ≈ 0.6 atol = 1.0e-12
    @test all(s -> insupport(d, s), rand(MersenneTwister(2), d, 20))
    # A direct keyword constructor is equivalent.
    d2 = TimeDelayedRedshiftDistribution(;
        params = _TD_PARAMS, source_frame_distribution = _td_ψ,
        time_delay_distribution = _TD_DELAY, maximum_redshift = 3.0, n_grid = 61
    )
    @test normalizer(d2) == normalizer(d)
end

@testset "time-delayed AD" begin
    # Derivative in the delay slope at α = -1 is finite.
    f = α -> time_delayed_merger_rate(
        1.0, _TD_PARAMS; formation_rate = _td_ψ, delay = PowerLawDelay(α, 0.02, 13.0)
    )
    g = ForwardDiff.derivative(f, -1.0)
    @test isfinite(g)
    @test g ≈ ForwardDiff.derivative(f, -1.0 + 1.0e-4) rtol = 1.0e-2

    h = Ωm -> time_delayed_merger_rate(
        1.0, merge(_TD_PARAMS, (Ωm = Ωm,)); formation_rate = _td_ψ, delay = _TD_DELAY
    )
    @test isfinite(ForwardDiff.derivative(h, 0.31))
end
