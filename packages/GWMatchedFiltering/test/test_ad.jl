@testset "Fixed-metric strain derivatives" begin
    x = [0.8, 0.2]
    for loss in (diagonal_loss, correlated_loss, statistics_loss)
        actual = ForwardDiff.gradient(loss, x)
        delta = 1.0e-5
        expected = [(loss(x + delta * e) - loss(x - delta * e)) / (2delta) for e in eachcol(Matrix{Float64}(I, 2, 2))]
        @test actual ≈ expected rtol = 1.0e-8 atol = 1.0e-8
    end
end

@testset "Plain-array whitening-factor strain derivatives" begin
    df = step(AD_F)
    W_diagonal = whitening_diagonal(AD_PSD)
    W_blocks = whitening_blocks(repeat(AD_S0, 1, 1, length(AD_F)))
    x = [0.8, 0.2]
    for W in (W_diagonal, W_blocks)
        loss(x) = log_likelihood_ratio(ad_signal(x), AD_DATA, W, df) + optimal_snr(ad_signal(x), W, df) +
            real(matched_filter_snr(ad_signal(x), AD_DATA, W, df)) + residual_squared(AD_DATA, ad_signal(x), W, df)
        actual = ForwardDiff.gradient(loss, x)
        delta = 1.0e-5
        expected = [(loss(x + delta * e) - loss(x - delta * e)) / (2delta) for e in eachcol(Matrix{Float64}(I, 2, 2))]
        @test actual ≈ expected rtol = 1.0e-8 atol = 1.0e-8
    end
end
