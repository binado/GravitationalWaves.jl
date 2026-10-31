using TestItems

@testsnippet Setup_whitened begin
end

@testitem "Plain-array inner products" setup = [Common, Setup_whitened] begin
    f = 0.0:0.5:2.0
    df = step(f)
    h = ComplexF64[0, 1, 2im, -1, 0]
    d = 2 .* h
    W = whitening_diagonal(reshape(fill(2.0, 5), :, 1))
    @test complex_inner_product(h, d, W, df) ≈ 12
    @test complex_inner_product(h, d, vec(W), df) ≈ 12
    @test complex_inner_product(reshape(h, :, 1), reshape(d, :, 1), W, df) ≈ 12
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

@testitem "Whitened contraction" setup = [Common, Setup_whitened] begin
    rng = Xoshiro(7)
    nf = 9
    h = randn(rng, ComplexF64, nf, 2)
    d = randn(rng, ComplexF64, nf, 2)
    Wd = whitening_diagonal(1 .+ rand(rng, nf, 2))
    @test GWMatchedFiltering.whitened_contraction(h, d, Wd) ≈ sum(conj.(h) .* Wd .^ 2 .* d)
    S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
    Wb = whitening_blocks(repeat(S0, 1, 1, nf))
    Wb[:, :, [1, nf]] .= 0
    # Oracle: explicit dense whitening per bin, then an ordinary dot.
    oracle = sum(dot(Wb[:, :, k] * h[k, :], Wb[:, :, k] * d[k, :]) for k in 1:nf)
    @test GWMatchedFiltering.whitened_contraction(h, d, Wb) ≈ oracle
    # Views take the same kernels as arrays.
    @test complex_inner_product(view(h, :, :), view(d, :, :), Wb, 0.25) ≈ complex_inner_product(h, d, Wb, 0.25)
    @test optimal_snr_squared(h .* 1.0e-22, Wd, 0.25) ≈ 1.0e-44 * optimal_snr_squared(h, Wd, 0.25)
    @test GWMatchedFiltering.filter_structure(Wd) isa GWMatchedFiltering.DiagonalFilter
    @test GWMatchedFiltering.filter_structure(Wb) isa GWMatchedFiltering.BlockFilter
end

@testitem "whiten" setup = [Common, Setup_whitened] begin
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
