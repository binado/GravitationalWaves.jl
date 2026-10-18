using FFTW

@testset "Physical Fourier convention" begin
    for T in (Float32, Float64), n in (2, 6, 16), nc in (1, 2)
        a = reshape(T.(sin.(1:(n * nc))), n, nc)
        nc == 1 && (a = vec(a))
        names = nc == 1 ? (:H1,) : (:H1, :L1)
        x = TimeSeries(a; epoch = 1126259460.125, sampling_frequency = 32, channels = names)
        h = to_frequencyseries(x)
        @test parent(h) ≈ rfft(a, 1) ./ 32
        @test epoch(h) == epoch(x)
        @test channels(h) == names
        @test nsamples(h) == n ÷ 2 + 1
        y = to_timeseries(h)
        @test parent(y) ≈ a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        @test check_same_grid(x, y) === nothing
        shifted = TimeSeries(a; epoch = 0, sampling_frequency = 32, channels = names)
        @test parent(to_frequencyseries(shifted)) == parent(h)
    end
    x = TimeSeries(collect(1.0:16); epoch = 0, sampling_frequency = 16, channels = (:H1,))
    h = to_frequencyseries(x)
    @test_throws ArgumentError to_frequencyseries(crop(x; start = 0, stop = 15 / 16))
    @test_throws ArgumentError to_frequencyseries(crop(x; start = 0, stop = 0))
    low_removed = crop(h; start = 1, stop = 9)
    @test_throws ArgumentError to_timeseries(low_removed)
    restored = pad(low_removed; before = 1)
    @test nsamples(to_timeseries(restored)) == 16
    high_removed = crop(h; start = 0, stop = 5)
    @test sampling_frequency(to_timeseries(high_removed)) == 8
    @test duration(to_timeseries(high_removed)) == 1
    @test sampling_frequency(to_timeseries(pad(high_removed; after = 4))) == 16
    @test_throws ArgumentError to_timeseries(crop(h; start = 0, stop = 1))
    dirty = copy(h)
    parent(dirty)[1] += 7im
    parent(dirty)[end] -= 9im
    @test parent(to_timeseries(dirty)) == parent(to_timeseries(h))
    @test parent(to_frequencyseries(TimeSeries(view(parent(x), :); epoch = 0, sampling_frequency = 16, channels = (:H1,)))) ≈ parent(h)
end
