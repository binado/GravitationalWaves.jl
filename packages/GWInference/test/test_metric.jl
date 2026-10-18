@testset "Setup, metadata, and frequency selection" begin
    f = 0.0:0.25:2.0
    names = (:H1, :L1, :V1)
    spectrum = fill(2.0, length(f), 3)
    metric = FrequencyDomainMetric(f, names, spectrum; minimum_frequency = 0.3, maximum_frequency = 1.6)
    @test metric.active_bins == 3:7
    @test metric.normalization == 1.0
    @test metric.metric isa DiagonalMetric
    spectrum .= 7.0
    @test all(==(0.5), metric.metric.weights)

    # Effective GWNoise PSDs have zeros at excluded bins; only the selected band is validated.
    spectrum .= 0
    spectrum[3:7, :] .= 2
    @test FrequencyDomainMetric(f, names, spectrum; minimum_frequency = 0.3, maximum_frequency = 1.6).active_bins == 3:7
    @test_throws ArgumentError FrequencyDomainMetric(f, names, spectrum)
    @test FrequencyDomainMetric(f, (:H1,), ones(9); nsamples = 16).active_bins == 2:8
    @test FrequencyDomainMetric(f, (:H1,), ones(9); nsamples = 17).active_bins == 2:9
    @test FrequencyDomainMetric(f, (:H1,), ones(9)).active_bins == 2:9
    @test_throws ArgumentError FrequencyDomainMetric(f, (:H1,), ones(9); nsamples = 14)
    @test_throws ArgumentError FrequencyDomainMetric(f, (:H1,), ones(9); minimum_frequency = 0, maximum_frequency = 0)
    @test_throws ArgumentError FrequencyDomainMetric(f, (:H1,), ones(9); minimum_frequency = -1)
    @test_throws ArgumentError FrequencyDomainMetric(f, (:H1,), ones(9); maximum_frequency = 3)
    @test_throws ArgumentError FrequencyDomainMetric(f, (:H1, "H1"), ones(9, 2))
    @test_throws ArgumentError FrequencyDomainStrain(f, (), zeros(9, 0))
    @test_throws DimensionMismatch FrequencyDomainStrain(f, names, zeros(3, 9))
    @test_throws DimensionMismatch FrequencyDomainMetric(f, names, ones(9))
    for invalid in (0.0, -1.0, Inf, NaN)
        s = ones(9)
        s[4] = invalid
        @test_throws ArgumentError FrequencyDomainMetric(f, (:H1,), s)
    end
    @test_throws ArgumentError FrequencyDomainStrain(2.0:-0.25:0.0, (:H1,), zeros(9))
    @test_throws ArgumentError FrequencyDomainStrain(range(0; length = 1), (:H1,), zeros(1))

    a = FrequencyDomainStrain(f, ("H1", "L1", "V1"), ones(9, 3))
    @test a.channels == names
    @test optimal_snr_squared(a, metric) == 7.5
    @test_throws ArgumentError optimal_snr(FrequencyDomainStrain(f, reverse(names), a.values), metric)
    @test_throws ArgumentError optimal_snr(FrequencyDomainStrain(0.0:0.5:4.0, names, a.values), metric)
    @test_throws DimensionMismatch optimal_snr(ones(8, 3), metric)
end

@testset "Correlated blocks and dense metric oracle" begin
    rng = Xoshiro(72)
    f = 0.0:0.5:4.0
    nc = 3
    S = zeros(ComplexF64, length(f), nc, nc)
    for k in 2:length(f)
        A = randn(rng, ComplexF64, nc, nc)
        S[k, :, :] = A * A' + I
    end
    metric = FrequencyDomainMetric(f, (:ET1, :ET2, :ET3), S)
    @test metric.metric isa BlockDiagonalMetric
    a = randn(rng, ComplexF64, length(f), nc)
    b = randn(rng, ComplexF64, length(f), nc)
    # Independent reference assembled in frequency-major flattening order.
    G = zeros(ComplexF64, nc * (length(f) - 1), nc * (length(f) - 1))
    for k in 2:length(f)
        indices = ((k - 2) * nc + 1):((k - 1) * nc)
        G[indices, indices] = 4step(f) * inv(S[k, :, :])
    end
    expected = dot(vec(permutedims(a[2:end, :])), G * vec(permutedims(b[2:end, :])))
    @test complex_inner_product(a, b, metric) ≈ expected rtol = 2.0e-14
    # Exercise the array-level kernel used by traced inputs as well as the host specialization.
    fallback = invoke(GWInference.contraction, Tuple{Any, Any, BlockDiagonalMetric}, a[2:end, :], b[2:end, :], metric.metric)
    @test metric.normalization * fallback ≈ expected rtol = 2.0e-14
    @test complex_inner_product(a, b, metric) ≈ conj(complex_inner_product(b, a, metric))
    @test optimal_snr_squared(a, metric) > 0
    @test abs(complex_inner_product(a, b, metric))^2 <= optimal_snr_squared(a, metric) * optimal_snr_squared(b, metric)

    order = [3, 1, 2]
    permuted = FrequencyDomainMetric(f, (:ET3, :ET1, :ET2), S[:, order, order])
    @test complex_inner_product(a[:, order], b[:, order], permuted) ≈ expected
    bad = copy(S)
    bad[2, 1, 2] += 1im
    @test_throws ArgumentError FrequencyDomainMetric(f, metric.channels, bad)
    bad[2, :, :] .= 0
    @test_throws PosDefException FrequencyDomainMetric(f, metric.channels, bad)
    bad[2, :, :] = Diagonal([1.0, 1.0, -1.0])
    @test_throws PosDefException FrequencyDomainMetric(f, metric.channels, bad)
    bad[2, 1, 1] = NaN
    @test_throws ArgumentError FrequencyDomainMetric(f, metric.channels, bad)
    @test_throws DimensionMismatch FrequencyDomainMetric(f, metric.channels, S[:, :, 1:2])

    psds = 1 .+ rand(rng, length(f), nc)
    diagonal_blocks = zeros(length(f), nc, nc)
    for i in 1:nc
        diagonal_blocks[:, i, i] = psds[:, i]
    end
    diagonal = FrequencyDomainMetric(f, metric.channels, psds)
    blocks = FrequencyDomainMetric(f, metric.channels, diagonal_blocks)
    @test complex_inner_product(a, b, diagonal) ≈ complex_inner_product(a, b, blocks)
end
