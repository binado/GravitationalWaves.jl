@testset "Series statistics and epoch alignment" begin
    f = 0:0.5:4
    names = (:H1, :L1)
    a = reshape(complex.(1:18), 9, 2)
    b = a .* (0.3 + 0.2im)
    for spectrum in (ones(9, 2), repeat(reshape([2.0 0.2; 0.2 1.0], 1, 2, 2), 9, 1, 1)), storage in (identity, x -> view(x, :, :))
        metric = FrequencyDomainMetric(f, names, spectrum)
        aa, bb = storage(a), storage(b)
        h = FrequencySeries(aa; frequencies = f, epoch = 1126259460.0, channels = names)
        d = FrequencySeries(bb; frequencies = f, epoch = epoch(h), channels = names)
        later = FrequencySeries(bb; frequencies = f, epoch = epoch(h) + 0.01, channels = names)
        reordered = FrequencySeries(bb; frequencies = f, epoch = epoch(h), channels = reverse(names))
        for statistic in (complex_inner_product, noise_weighted_inner_product, matched_filter_snr, residual_squared, likelihood_statistics, log_likelihood_ratio)
            @test statistic(h, d, metric) == statistic(aa, bb, metric)
            @test_throws ArgumentError statistic(h, later, metric)
            @test_throws ArgumentError statistic(h, reordered, metric)
            @test_throws ArgumentError statistic(h, bb, metric)
            @test_throws ArgumentError statistic(aa, d, metric)
        end
        @test optimal_snr(h, metric) == optimal_snr(aa, metric)
        @test optimal_snr_squared(later, metric) == optimal_snr_squared(bb, metric)
        @test GWMatchedFiltering.check_array(h, metric) === aa
        @test matched_filter_snr(h, d, metric; check = false) == matched_filter_snr(aa, bb, metric; check = false)
    end
    metric = FrequencyDomainMetric(f, (:H1,), ones(9))
    h = FrequencySeries(vec(a[:, 1]); frequencies = f, epoch = 0, channels = (:H1,))
    @test optimal_snr(h, metric) == optimal_snr(parent(h), metric)
    @test ForwardDiff.derivative(scale -> optimal_snr_squared(with_data(h, scale .* parent(h)), metric), 1.0) ≈ 2optimal_snr_squared(h, metric)
end
