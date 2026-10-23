include("whitening_helpers.jl")

@testset "Plain-array inner products" begin
    f = 0.0:0.5:2.0
    df = step(f)
    h = ComplexF64[0, 1, 2im, -1, 0]
    d = 2 .* h
    W = whitening_diagonal(reshape(fill(2.0, 5), :, 1))
    @test complex_inner_product(h, d, W, df) ≈ 12
    @test complex_inner_product(h, d, vec(W), df) ≈ 12
    @test noise_weighted_inner_product(h, d, W, df) ≈ 12
    @test optimal_snr_squared(h, W, df) ≈ 6
    @test optimal_snr(h, W, df) ≈ sqrt(6)
    @test matched_filter_snr(h, d, W, df) ≈ 2sqrt(6)
    @test matched_filter_snr(h, im .* h, W, df) ≈ im * sqrt(6)
    @test optimal_snr(zero(h), W, df) == 0
    @test_throws DomainError check_template(zero(h), W, df)
    @test_throws DimensionMismatch check_template(ComplexF64[1 1; 1 1; 1 1; 1 1; 1 1], W, df)
    @test check_template(h, W, df) === nothing
    @test isnan(matched_filter_snr(zero(h), d, W, df))
    @test residual_squared(d, h, W, df) ≈ 6
    stats = likelihood_statistics(h, d, W, df)
    @test stats.h_inner_d ≈ 12 && stats.h_inner_h ≈ 6
    @test log_likelihood_ratio(h, d, W, df) ≈ 9
    @test log_likelihood_ratio(h, d, W, df) ≈ (optimal_snr_squared(d, W, df) - residual_squared(d, h, W, df)) / 2
    # Mixed vector/matrix single-channel arguments.
    @test complex_inner_product(h, reshape(d, :, 1), W, df) ≈ 12
    @test complex_inner_product(reshape(h, :, 1), d, W, df) ≈ 12

    # Zero weights mark inactive bins.
    W0 = copy(W)
    W0[1:2] .= 0
    @test optimal_snr_squared(h, W0, df) ≈ 4df * (4 / 2 + 1 / 2)

    # Equicorrelated three-channel noise: the common-mode eigenvalue is 1+2r.
    r = 0.2
    S = repeat(reshape(fill(r, 3, 3) + (1 - r) * I, 3, 3, 1), 1, 1, length(f))
    Wb = whitening_blocks(S)
    common = ones(ComplexF64, length(f), 3)
    @test optimal_snr_squared(common, Wb, df) ≈ 4df * 5 * 3 / (1 + 2r)
    differential = repeat(reshape([1.0, -1.0, 0.0], 1, 3), length(f), 1)
    @test optimal_snr_squared(differential, Wb, df) ≈ 4df * 5 * 2 / (1 - r)
    # One channel in (1, 1, nf) block storage matches the diagonal path.
    Wone = reshape(vec(W), 1, 1, :)
    @test complex_inner_product(h, d, Wone, df) ≈ 12
    @test complex_inner_product(reshape(h, :, 1), reshape(d, :, 1), Wone, df) ≈ 12
end

@testset "Plain arrays agree with the metric API" begin
    rng = Xoshiro(7)
    f = 0.0:0.25:2.0
    nf = length(f)
    h = randn(rng, ComplexF64, nf, 2)
    d = randn(rng, ComplexF64, nf, 2)
    psd = 1 .+ rand(rng, nf, 2)
    metric_domain = frequency_prototype(; frequencies = f, epoch = 0.0, channels = (:A, :B))
    diagonal = FrequencyDomainMetric(metric_domain, psd; minimum_frequency = f[2], maximum_frequency = f[end - 1])
    Wd = whitening_diagonal(psd)
    Wd[[1, nf], :] .= 0
    @test complex_inner_product(h, d, Wd, step(f)) ≈ complex_inner_product(h, d, diagonal) rtol = 1.0e-13

    S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
    S = repeat(reshape(S0, 1, 2, 2), nf, 1, 1)
    blocks = FrequencyDomainMetric(metric_domain, S; minimum_frequency = f[2], maximum_frequency = f[end - 1])
    Wb = whitening_blocks(permutedims(S, (2, 3, 1)))
    Wb[:, :, [1, nf]] .= 0
    @test complex_inner_product(h, d, Wb, step(f)) ≈ complex_inner_product(h, d, blocks) rtol = 1.0e-13
    @test matched_filter_snr(h, d, Wb, step(f)) ≈ matched_filter_snr(h, d, blocks) rtol = 1.0e-13
    # Views take the same kernels as arrays; the broadcast fallback agrees with them.
    @test complex_inner_product(view(h, :, :), view(d, :, :), Wb, step(f)) ≈ complex_inner_product(h, d, Wb, step(f))
    @test GWMatchedFiltering.whitened_contraction(h, d, @view(Wb[:, :, :])) ≈ GWMatchedFiltering.whitened_contraction(h, d, Wb)
    @test GWMatchedFiltering.whitened_contraction(h, d, Wd) ≈ sum(conj.(h) .* Wd .^ 2 .* d)
    @test @allocated(complex_inner_product(h, d, Wd, 0.25)) == 0
    @test @allocated(complex_inner_product(h, d, Wb, 0.25)) == 0
end

@testset "whiten" begin
    rng = Xoshiro(3)
    nf = 6
    d = randn(rng, ComplexF64, nf, 2)
    Wd = rand(rng, nf, 2)
    @test whiten(d, Wd) == Wd .* d
    @test whiten(d[:, 1], Wd[:, 1]) == Wd[:, 1] .* d[:, 1]
    @test whiten(d[:, 1], Wd[:, 1:1]) == Wd[:, 1] .* d[:, 1]
    S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
    Wb = whitening_blocks(repeat(S0, 1, 1, nf))
    white = whiten(d, Wb)
    @test white[3, :] ≈ Wb[:, :, 3] * d[3, :]
    out = similar(d)
    @test whiten!(out, d, Wb) === out && out ≈ white
    # Whitened Euclidean sums reproduce the inner product.
    @test 4 * 0.5 * sum(conj.(white) .* whiten(d, Wb)) ≈ complex_inner_product(d, d, Wb, 0.5)
end
