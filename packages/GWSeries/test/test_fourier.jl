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
        full = fft(x)
        @test metadata(full).spectrum == :twosided
        @test frequencies(full) == range(-fld(n, 2) / (n * 0.125); step = inv(n * 0.125), length = n)
        @test parent(full) ≈ fftshift(fft(a, 1), 1) .* 0.125
        @test parent(ifft(full)) ≈ a
        @test parent(to_timedomain(to_frequencydomain(x); nsamples = n)) ≈ a
        @test_throws ArgumentError rfft(x, 2)
        @test_throws ArgumentError ifft(h)
        @test_throws ArgumentError irfft(full, n)
        @test_throws ArgumentError fft!(x)
        @test_throws ArgumentError to_timedomain(h[2:end, :])
        for real_transform in (true, false)
            p = plan_to_frequencydomain(x; real_transform)
            actual = transform(p, x)
            expected = real_transform ? h : full
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
    @test parent(to_timedomain(to_frequencydomain(x))) ≈ parent(x)
    h = to_frequencydomain(real.(x))
    @test sum(abs2, parent(x)) * sampling_interval(x) ≈ sum(abs2, parent(fft(x))) * frequency_spacing(fft(x))
    @test_throws ArgumentError to_frequencydomain(TimeSeries(ones(1); epoch = 0, sampling_frequency = 1, channels = (:H1,)))
    @test_throws DimensionMismatch transform!(zeros(1, 1), plan_to_frequencydomain(real.(x)), parent(real.(x)))
    @test to_timedomain(h) isa TimeSeries
end
