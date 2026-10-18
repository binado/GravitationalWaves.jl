using GWSeries

@testset "Series noise delegation" begin
    curve = TabulatedPSD([0, 64], [2, 2])
    even = GaussianNoise(curve; sampling_frequency = 128, nsamples = 128)
    odd = GaussianNoise(curve; sampling_frequency = 128, nsamples = 127)
    h = noise_fd(Xoshiro(42), FrequencySeries, even; epoch = 1126259460.0, channels = (:H1,))
    @test parent(h) == noise_fd(Xoshiro(42), even)
    @test epoch(h) == 1126259460.0
    @test parent(to_timeseries(h)) ≈ noise_td(Xoshiro(42), even)
    @test noise_fd!(Xoshiro(43), h, even) === h
    @test parent(h) == noise_fd(Xoshiro(43), even)
    @test noise_fd!(h, even) === h
    @test noise_fd(FrequencySeries, even; epoch = 0, channels = (:H1,)) isa FrequencySeries
    for model in (even, odd)
        x = noise_td(Xoshiro(42), TimeSeries, model; epoch = 1126259460.0, channels = (:H1,))
        @test parent(x) == noise_td(Xoshiro(42), model)
        @test noise_td!(Xoshiro(43), x, model) === x
        @test parent(x) == noise_td(Xoshiro(43), model)
        @test noise_td!(x, model) === x
        @test noise_td(TimeSeries, model; epoch = 0, channels = (:H1,)) isa TimeSeries
    end
    @test_throws ArgumentError noise_fd(FrequencySeries, odd; epoch = 0, channels = (:H1,))
    @test_throws DimensionMismatch noise_fd(FrequencySeries, even; epoch = 0, channels = (:H1, :L1))
    bad = FrequencySeries(ones(ComplexF64, 65); frequencies = 1:1:65, epoch = 0, channels = (:H1,))
    saved = copy(parent(bad))
    rng = Xoshiro(41)
    control = copy(rng)
    @test_throws ArgumentError noise_fd!(rng, bad, even)
    @test parent(bad) == saved
    @test rand(rng) == rand(control)
    badtime = TimeSeries(ones(128); epoch = 0, sampling_frequency = 256, channels = (:H1,))
    @test_throws ArgumentError noise_td!(badtime, even)
    @test all(isone, parent(badtime))
    matrix = FrequencySeries(zeros(ComplexF64, 65, 1); frequencies = frequencies(even), epoch = 0, channels = (:H1,))
    @test_throws ArgumentError noise_fd!(matrix, even)
end
