using TestItems

@testsnippet Setup_arrays begin
    using LinearAlgebra
end

@testitem "Plain-array noise covariance" setup = [Common, Setup_arrays] begin
    @testset "band_psd" begin
        curve = TabulatedPSD([0, 64], [2, 2])
        f = range(0.0; step = 1.0, length = 65)
        S = band_psd(curve, f; fmin = 10.5, fmax = 30)
        @test S == [10.5 <= x <= 30 ? 2.0 : 0.0 for x in f]
        full = band_psd(curve, f)
        @test full[1] == 0 && full[end] == 0 && all(==(2.0), full[2:(end - 1)])
        @test_throws ArgumentError band_psd(curve, f; fmax = 65)
        @test_throws ArgumentError band_psd(curve, f; fmin = 40, fmax = 10)
    end

    @testset "whitening_factor" begin
        S = Float64[0 0; 2 8; 4 2; 0 0]
        W = whitening_factor(S)
        @test W == [0 0; 1 / sqrt(2) 1 / sqrt(8); 0.5 1 / sqrt(2); 0 0]
        @test_throws ArgumentError whitening_factor([1.0 -1.0])
        @test_throws ArgumentError whitening_factor([1.0 NaN])

        S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
        S3 = zeros(ComplexF64, 2, 2, 4)
        S3[:, :, 2] = S0
        S3[:, :, 3] = 3S0
        W3 = whitening_factor(S3)
        L = cholesky(Hermitian(S0)).L
        @test W3[:, :, 2] ≈ inv(L)
        @test W3[:, :, 3] ≈ inv(sqrt(3) * L)
        @test iszero(W3[:, :, 1]) && iszero(W3[:, :, 4])
        @test W3[1, 2, 2] == 0  # lower triangular
        for k in 2:3
            @test W3[:, :, k] * S3[:, :, k] * W3[:, :, k]' ≈ I atol = 1.0e-14
        end
        bad = copy(S3)
        bad[:, :, 2] = ComplexF64[1 1; 1 1]
        @test_throws ArgumentError whitening_factor(bad)
        bad[:, :, 2] = ComplexF64[1 1im; 1im 1]
        @test_throws ArgumentError whitening_factor(bad)
        bad[:, :, 2] = ComplexF64[1 NaN; NaN 1]
        @test_throws ArgumentError whitening_factor(bad)
        @test_throws DimensionMismatch whitening_factor(zeros(2, 3, 4))
        @test whitening_factor(reshape(S[:, 2], 1, 1, :)) ≈ reshape(whitening_factor(S[:, 2:2]), 1, 1, :)
    end

    @testset "sample_noise and sample_noise!" begin
        curve = TabulatedPSD([0, 64], [1, 129])
        S = reshape(band_psd(curve, range(0.0; step = 1.0, length = 65)), :, 1)
        df = 1.0
        draw = sample_noise(Xoshiro(5), S, df)
        @test size(draw) == size(S) && eltype(draw) == ComplexF64
        @test draw == sample_noise(Xoshiro(5), S, df)
        @test draw != sample_noise(Xoshiro(6), S, df)
        @test all(iszero, draw[[1, end], :])
        out = similar(draw)
        @test sample_noise!(Xoshiro(5), out, S, df) === out && out == draw
        @test size(sample_noise(S, df)) == size(S)
        @test_throws DimensionMismatch sample_noise!(zeros(ComplexF64, 3, 1), S, df)
        @test_throws ArgumentError sample_noise(Xoshiro(1), [1.0 2.0], 0.0)
    end

    @testset "Sampler covariance is S/(2df)" begin
        rng = Xoshiro(2026)
        df = 0.25
        draws = 40_000
        # Independent channels with different spectra.
        Sd = [0.0 0.0; 2.0 8.0; 4.5 0.5]
        acc = zeros(ComplexF64, size(Sd))
        for _ in 1:draws
            acc .+= abs2.(sample_noise(rng, Sd, df))
        end
        @test acc[1, :] == [0, 0]
        @test (acc ./ draws)[2:3, :] ≈ Sd[2:3, :] ./ (2df) rtol = 0.03
        # Correlated channels: E[n n†] = S/(2df) per bin, E[n nᵀ] = 0.
        S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
        Sb = zeros(ComplexF64, 2, 2, 3)
        Sb[:, :, 2] = S0
        Sb[:, :, 3] = 2S0
        cov = zeros(ComplexF64, 2, 2, 3)
        pseudo = zeros(ComplexF64, 2, 2, 3)
        for _ in 1:draws
            n = sample_noise(rng, Sb, df)
            for k in 1:3
                cov[:, :, k] .+= n[k, :] * n[k, :]'
                pseudo[:, :, k] .+= n[k, :] * transpose(n[k, :])
            end
        end
        cov ./= draws
        pseudo ./= draws
        @test iszero(cov[:, :, 1])
        @test cov[:, :, 2] ≈ S0 / (2df) rtol = 0.04
        @test cov[:, :, 3] ≈ 2S0 / (2df) rtol = 0.04
        @test maximum(abs, pseudo[:, :, 2]) < 0.05 * maximum(abs, S0 / (2df))
        @test_throws DimensionMismatch sample_noise!(zeros(ComplexF64, 2, 2), Sb, df)
    end
end

@testitem "Noise and matched-filtering normalization contract" setup = [Common, Setup_arrays] begin
    using GWMatchedFiltering: whiten, complex_inner_product, noise_weighted_inner_product, optimal_snr_squared
    rng = Xoshiro(11)
    df = 0.25
    draws = 20_000
    # Whitening drawn noise gives E[(W n)(W n)†] = I/(2df) per active bin, independent of the PSD.
    Sd = [0.0 0.0; 2.0 8.0; 4.5 0.5]
    Wd = whitening_factor(Sd)
    S0 = ComplexF64[2 0.3 + 0.4im; 0.3 - 0.4im 1]
    Sb = zeros(ComplexF64, 2, 2, 3)
    Sb[:, :, 2] = S0
    Sb[:, :, 3] = 2S0
    Wb = whitening_factor(Sb)
    for (S, W, nactive) in ((Sd, Wd, 4), (Sb, Wb, 4))
        energy = 0.0
        for _ in 1:draws
            n = sample_noise(rng, S, df)
            white = whiten(n, W)
            @test size(white) == size(n)
            energy += sum(abs2, white)
            # Euclidean sum of whitened components reproduces the 4df-weighted product.
        end
        # Each active (bin, channel) component carries E|z|² = 1/(2df).
        @test energy / draws ≈ nactive / (2df) rtol = 0.01
    end
    # ⟨n, n⟩ = 4df Σ |W n|² has mean 2 per active complex component (chi-squared, 2 dof).
    for (S, W, nactive) in ((Sd, Wd, 4), (Sb, Wb, 4))
        mean_snr2 = sum(optimal_snr_squared(sample_noise(rng, S, df), W, df) for _ in 1:draws) / draws
        @test mean_snr2 ≈ 2nactive rtol = 0.01
    end
    # The matched-filter statistic ⟨h, n⟩/√⟨h, h⟩ has E|ρ|² = 2 (unit variance in Re and Im).
    h = ones(ComplexF64, 3, 2)
    norm = sqrt(optimal_snr_squared(h, Wb, df))
    rho = [complex_inner_product(h, sample_noise(rng, Sb, df), Wb, df) / norm for _ in 1:draws]
    @test sum(abs2, rho) / draws ≈ 2 rtol = 0.03
    @test abs(sum(rho) / draws) < 0.05
end
