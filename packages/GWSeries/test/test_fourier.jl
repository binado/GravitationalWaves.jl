@testset "Physical AbstractFFTs transforms" begin
    rng = MersenneTwister(37)
    for T in (Float32, Float64), n in (7, 8, 17, 32), nc in (1, 3)
        a = randn(rng, T, n, nc)
        x = TimeSeries(a; epoch = 1126259460.125, sampling_interval = 0.125, channels = Tuple(Symbol("c", i) for i in 1:nc))
        h = rfft(x)
        @test h isa FrequencySeries
        @test parent(h) ≈ rfft(a, 1) .* 0.125
        @test frequency_spacing(h) == inv(n * 0.125)
        @test epoch(h) == epoch(x)
        @test channels(h) == channels(x)
        @test parent(irfft(h, n)) ≈ a
        @test metadata(irfft(h, n)) == metadata(x)
        @test parent(to_timedomain(to_frequencydomain(x); nsamples = n)) ≈ a
        @test_throws ArgumentError rfft(x, 2)
        @test_throws ArgumentError ifft(h)
        @test_throws ArgumentError fft(x)
        @test_throws ArgumentError fft!(x)
        @test_throws MethodError to_timedomain(h[2:end, :]; nsamples = n)
        @test_throws UndefKeywordError to_timedomain(h)
        @test_throws UndefKeywordError plan_to_timedomain(h)
        @test_throws ArgumentError irfft(h, 2n)
        let
            p = plan_to_frequencydomain(x)
            actual = transform(p, x)
            expected = h
            @test parent(actual) ≈ parent(expected)
            @test metadata(p) == metadata(expected)
            raw = transform(p, parent(x))
            @test raw ≈ parent(expected)
            target = copy(expected)
            transform!(target, p, x)
            @test parent(target) ≈ parent(expected)
            @test transform!(parent(target), p, x) === parent(target)
            @test transform!(target, p, parent(x)) === target
            wrong_epoch = FrequencySeries(copy(parent(target)); metadata(target)..., epoch = epoch(target) + 1)
            @test_throws ArgumentError transform!(wrong_epoch, p, parent(x))
            shifted_input = TimeSeries(copy(parent(x)); metadata(x)..., epoch = epoch(x) + 1)
            @test_throws ArgumentError transform!(parent(target), p, shifted_input)
            @test_throws ArgumentError transform(p, expected)
            inverse = plan_to_timedomain(actual; nsamples = n)
            saved = copy(parent(actual))
            @test parent(transform(inverse, actual)) ≈ a
            @test parent(actual) == saved
            @test transform(inverse, parent(actual)) ≈ a
        end
    end
    x = TimeSeries(randn(rng, 16, 2) .+ im .* randn(rng, 16, 2); epoch = -2, sampling_frequency = 8, channels = (:plus, :cross))
    @test_throws ArgumentError rfft(x)
    @test_throws ArgumentError plan_fft(x)
    @test_throws ArgumentError to_frequencydomain(x)
    @test_throws ArgumentError plan_to_frequencydomain(x)
    h = to_frequencydomain(real.(x))
    @test sum(abs2, parent(x)) ≈ sum(abs2, fft(parent(x), 1)) / size(x, 1)
    @test_throws ArgumentError to_frequencydomain(TimeSeries(ones(1); epoch = 0, sampling_frequency = 1, channels = (:H1,)))
    @test_throws DimensionMismatch transform!(zeros(1, 1), plan_to_frequencydomain(real.(x)), parent(real.(x)))
    @test to_timedomain(h; nsamples = 16) isa TimeSeries
    bins = frequency_bins(h; fmin = 1.1, fmax = 2.6)
    full = copy(h)
    fill!(full, 0)
    full[bins, :] = h[bins, :]
    reference = zeros(eltype(h), size(h))
    reference[bins, :] = parent(h)[bins, :]
    @test parent(irfft(full, 16)) ≈ irfft(reference, 16, 1) .* (16 * frequency_spacing(h))
end
