@testset "Explicit domain validation and unchecked evaluation" begin
    f = 0:0.5:4
    names = (:H1, :L1)
    fd = FrequencyDomain(; frequencies = f, epoch = 1126259460.0, channels = names)
    later = FrequencyDomain(; frequencies = f, epoch = 1126259461.0, channels = names)
    reordered = FrequencyDomain(; frequencies = f, epoch = epoch(fd), channels = reverse(names))
    wider = FrequencyDomain(; frequencies = 0:1:8, epoch = epoch(fd), channels = names)
    longer = FrequencyDomain(; frequencies = 0:0.5:4.5, epoch = epoch(fd), channels = names)
    a = reshape(complex.(1:18), 9, 2);b = a .* (0.3 + 0.2im)
    for spectrum in (ones(9, 2), repeat(reshape([2.0 0.2;0.2 1.0], 1, 2, 2), 9, 1, 1)), storage in (identity, x -> view(x, :, :))
        metric = FrequencyDomainMetric(fd, spectrum)
        aa, bb = storage(a), storage(b)
        @test check_data(fd, aa) === nothing
        @test check_compatible(fd, metric) === nothing
        @test check_compatible(later, metric) === nothing
        @test_throws ArgumentError check_same_grid(fd, later)
        @test_throws ArgumentError check_compatible(reordered, metric)
        @test_throws ArgumentError check_compatible(wider, metric)
        @test check_template(fd, aa, metric) === nothing
        @test_throws DomainError check_template(fd, zero(aa), metric)
        @test_throws DomainError check_template(fd, fill(ComplexF64(NaN), size(aa)), metric)
        for statistic in (complex_inner_product, noise_weighted_inner_product, matched_filter_snr, residual_squared, likelihood_statistics, log_likelihood_ratio)
            @test statistic(fd, aa, bb, metric) == statistic(later, aa, bb, metric)
            # Deliberately incompatible labels would throw if a validator ran here.
            @test statistic(reordered, aa, bb, metric) == statistic(fd, aa, bb, metric)
            @test statistic(longer, aa, bb, metric) == statistic(fd, aa, bb, metric)
        end
        @test complex_inner_product(wider, aa, bb, metric) == 2complex_inner_product(fd, aa, bb, metric)
        @test ForwardDiff.derivative(scale -> optimal_snr_squared(fd, scale .* aa, metric), 1.0) ≈ 2optimal_snr_squared(fd, aa, metric)
    end
end

@testset "Single-channel vectors, views, and column matrices" begin
    fd = FrequencyDomain(; frequencies = 0:0.5:4, epoch = 0, channels = (:H1,))
    backing = complex.(sin.(1:18), cos.(1:18))
    for a in (backing[1:9], view(backing, 1:2:18)), spectrum in (fill(2.0, 9), fill(2.0 + 0im, 9, 1, 1))
        b = (0.4 + 0.3im) .* a
        ma, mb = reshape(collect(a), :, 1), reshape(b, :, 1)
        metric = FrequencyDomainMetric(fd, spectrum; nsamples = 16)
        expected = complex_inner_product(fd, ma, mb, metric)
        @test complex_inner_product(fd, a, b, metric) ≈ expected
        @test complex_inner_product(fd, a, mb, metric) ≈ expected
        @test complex_inner_product(fd, ma, b, metric) ≈ expected
        @test matched_filter_snr(fd, a, b, metric) ≈ matched_filter_snr(fd, ma, mb, metric)
        @test residual_squared(fd, a, b, metric) ≈ residual_squared(fd, ma, mb, metric)
    end
end
