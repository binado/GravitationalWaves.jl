if !isdefined(@__MODULE__, :time_prototype)
    include(joinpath(@__DIR__, "../../GWSeries/test/prototypes.jl"))
end
@testset "Checked series and unchecked array statistics" begin
    f = 0:0.5:4
    names = (:H1, :L1)
    a = FrequencySeries(reshape(complex.(1:18), 9, 2); frequency_spacing = step(f), epoch = 1126259460.0, channels = names)
    b = (0.3 + 0.2im) .* a
    metric = FrequencyDomainMetric(a, ones(9, 2))
    bins = frequency_bins(a; fmin = 1.1, fmax = 3.1)
    band_metric = FrequencyDomainMetric(ones(length(bins), 2); frequencies = frequencies(a)[bins], channels = names)
    full_metric = FrequencyDomainMetric(a, ones(9, 2); minimum_frequency = 1.1, maximum_frequency = 3.1)
    @test complex_inner_product(a[bins, :], b[bins, :], band_metric) ≈ complex_inner_product(a, b, full_metric)
    @test check_template(a, metric) === nothing
    @test check_template(parent(a), metric) === nothing
    @test_throws DomainError check_template(zero(parent(a)), metric)
    later = FrequencySeries(copy(parent(a)); metadata(a)..., epoch = epoch(a) + 1)
    reordered = FrequencySeries(copy(parent(a)); metadata(a)..., channels = reverse(names))
    @test check_compatible(later, metric) === nothing
    for statistic in (complex_inner_product, noise_weighted_inner_product, matched_filter_snr, residual_squared, likelihood_statistics, log_likelihood_ratio)
        @test statistic(a, b, metric) == statistic(parent(a), parent(b), metric)
        @test_throws ArgumentError statistic(a, later, metric)
        @test_throws ArgumentError statistic(a, reordered, metric)
    end

end

@testset "Single-channel vectors, views, and column matrices" begin
    fd = frequency_prototype(; frequencies = 0:0.5:4, epoch = 0, channels = (:H1,))
    backing = complex.(sin.(1:18), cos.(1:18))
    for a in (backing[1:9], view(backing, 1:2:18)), spectrum in (fill(2.0, 9), fill(2.0 + 0im, 9, 1, 1))
        b = (0.4 + 0.3im) .* a
        ma, mb = reshape(collect(a), :, 1), reshape(b, :, 1)
        metric = FrequencyDomainMetric(fd, spectrum; nsamples = 16)
        expected = complex_inner_product(ma, mb, metric)
        @test complex_inner_product(a, b, metric) ≈ expected
        @test complex_inner_product(a, mb, metric) ≈ expected
        @test complex_inner_product(ma, b, metric) ≈ expected
        @test matched_filter_snr(a, b, metric) ≈ matched_filter_snr(ma, mb, metric)
        @test residual_squared(a, b, metric) ≈ residual_squared(ma, mb, metric)
    end
end
