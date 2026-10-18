@testset "SNRs and Gaussian statistics" begin
    f = 0.0:0.5:2.0
    metric = FrequencyDomainMetric(f, (:H1,), fill(2.0, 5))
    h = ComplexF64[0, 1, 2im, -1, 0]
    d = 2 .* h
    @test complex_inner_product(h, d, metric) == 12
    @test noise_weighted_inner_product(h, d, metric) == 12
    @test optimal_snr_squared(h, metric) == 6
    @test optimal_snr(h, metric) == sqrt(6)
    @test matched_filter_snr(h, d, metric) ≈ 2sqrt(6)
    @test matched_filter_snr(h, d, metric; check = false) ≈ 2sqrt(6)
    @test matched_filter_snr(h, im .* h, metric) ≈ im * sqrt(6)
    @test optimal_snr(zero(h), metric) == 0
    @test_throws DomainError matched_filter_snr(zero(h), d, metric)
    @test residual_squared(d, h, metric) == 6
    @test likelihood_statistics(h, d, metric) == (; h_inner_d = 12 + 0im, h_inner_h = 6)
    @test log_likelihood_ratio(h, d, metric) == 9
    @test log_likelihood_ratio(h, d, metric) ≈ (optimal_snr_squared(d, metric) - residual_squared(d, h, metric)) / 2
    @test optimal_snr_squared(h .* 1.0e-22, metric) ≈ 6.0e-44

    # Three channels with equicorrelated noise: the common-mode eigenvalue is 1+2r.
    r = 0.2
    S = repeat(reshape(fill(r, 3, 3) + (1 - r) * I, 1, 3, 3), length(f), 1, 1)
    correlated = FrequencyDomainMetric(f, (:ET1, :ET2, :ET3), S)
    common = ones(ComplexF64, length(f), 3)
    @test optimal_snr_squared(common, correlated) ≈ 4step(f) * 4 * 3 / (1 + 2r)
    differential = repeat(reshape([1.0, -1.0, 0.0], 1, 3), length(f), 1)
    @test optimal_snr_squared(differential, correlated) ≈ 4step(f) * 4 * 2 / (1 - r)
end

@testset "Correlated Gaussian noise normalization" begin
    rng = Xoshiro(2026)
    f = 0.0:0.25:2.0
    S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
    S = repeat(reshape(S0, 1, 2, 2), length(f), 1, 1)
    metric = FrequencyDomainMetric(f, (:A, :B), S; nsamples = 16)
    L = cholesky(Hermitian(S0)).L
    h = ones(ComplexF64, length(f), 2)
    rho = Vector{ComplexF64}(undef, 4000)
    energy = 0.0
    for draw in eachindex(rho)
        # E[n_k n_k'] = S_k/(2Δf), matching GWNoise's physical Fourier units.
        n = permutedims(L * randn(rng, ComplexF64, 2, length(f))) / sqrt(2step(f))
        rho[draw] = matched_filter_snr(h, n, metric)
        energy += optimal_snr_squared(n, metric)
    end
    @test abs(sum(rho) / length(rho)) < 0.08
    @test sum(abs2, rho) / length(rho) ≈ 2 rtol = 0.06
    @test energy / length(rho) ≈ 2 * 2 * length(metric.active_bins) rtol = 0.03
end
