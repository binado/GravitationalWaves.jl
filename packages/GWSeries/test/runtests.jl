using Test
using GWSeries

@testset "Storage and metadata" begin
    a = Float32[1, 2, 3, 4]
    x = TimeSeries(a; epoch = 1126259460.125, sampling_frequency = 4, channels = ("H1",))
    @test parent(x) === a
    @test epoch(x) === 1126259460.125
    @test channels(x) == (:H1,)
    @test sample_spacing(x) == 0.25
    @test sampling_frequency(x) == 4
    @test duration(x) == 1
    @test times(x) == range(epoch(x); step = 0.25, length = 4)
    @test !(x isa AbstractArray)
    a[1] = 5
    @test parent(x)[1] == 5
    y = copy(x)
    @test parent(y) !== a
    @test parent(y) == a
    b = fill(1im, 4)
    @test parent(with_data(x, b)) === b
    @test_throws DimensionMismatch with_data(x, ones(4, 1))
    @test_throws DimensionMismatch TimeSeries(a; epoch = 0, sampling_frequency = 4, channels = (:H1, :L1))
    for names in ((), (:H1, "H1"), (1,))
        @test_throws ArgumentError TimeSeries(ones(4, 2); epoch = 0, sampling_frequency = 4, channels = names)
    end
    for value in (0, -1, NaN, Inf)
        @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_frequency = value, channels = (:H1,))
    end
    @test_throws ArgumentError TimeSeries(a; epoch = Inf, sampling_frequency = 4, channels = (:H1,))
    big = TimeSeries(a; epoch = big"1126259460.123456789", sampling_frequency = big"4", channels = (:H1,))
    @test epoch(big) isa BigFloat
    @test sample_spacing(big) isa BigFloat
    @test check_same_grid(big, copy(big)) === nothing
    h = FrequencySeries(a; frequencies = 2:0.25:2.75, epoch = epoch(x), channels = (:H1,))
    @test parent(h) === a
    @test frequencies(h) == 2:0.25:2.75
    @test frequency_spacing(h) == 0.25
    @test nsamples(h) == 4
    @test_throws MethodError duration(h)
    @test_throws ArgumentError check_same_grid(x, h)
    @test_throws DimensionMismatch FrequencySeries(a; frequencies = 0:1:2, epoch = 0, channels = (:H1,))
    for f in (-1:1:2, 4:-1:1, range(0.0; step = Inf, length = 4))
        @test_throws ArgumentError FrequencySeries(a; frequencies = f, epoch = 0, channels = (:H1,))
    end
    singleton = FrequencySeries([1]; frequencies = range(3; step = 2, length = 1), epoch = 0, channels = (:H1,))
    @test frequency_spacing(singleton) == 2
end

@testset "Grid operations and compatibility" begin
    t0 = 1126259460.1
    for a in (collect(1.0:8), reshape(collect(1.0:16), 8, 2))
        names = ndims(a) == 1 ? (:H1,) : (:H1, :L1)
        x = TimeSeries(a; epoch = t0, sampling_frequency = 10, channels = names)
        y = crop(x; start = t0 + 0.2, stop = t0 + 0.6)
        @test nsamples(y) == 4
        @test epoch(y) == t0 + 0.2
        @test parent(y) == (ndims(a) == 1 ? a[3:6] : a[3:6, :])
        @test parent(y) !== parent(x)
        @test check_same_grid(x, copy(x)) === nothing
        @test_throws ArgumentError crop(x; start = t0 + 0.21, stop = t0 + 0.6)
        @test_throws ArgumentError crop(x; start = t0 - 0.1, stop = t0)
        @test_throws ArgumentError crop(x; start = t0 + 0.6, stop = t0 + 0.2)
        empty = crop(x; start = t0 + 0.8, stop = t0 + 0.8)
        @test nsamples(empty) == 0
        @test epoch(empty) == t0 + 0.8
        padded = pad(x; before = 2, after = 3)
        @test nsamples(padded) == 13
        @test epoch(padded) == t0 - 0.2
        @test check_same_grid(crop(padded; start = t0, stop = t0 + 0.8), x) === nothing
        @test_throws ArgumentError pad(x; before = -1)
        h = FrequencySeries(a; frequencies = 2:0.5:5.5, epoch = t0, channels = names)
        hp = pad(h; before = 4, after = 2)
        @test first(frequencies(hp)) == 0
        @test epoch(hp) == t0
        @test parent(crop(hp; start = 2, stop = 6)) == a
        @test_throws ArgumentError pad(h; before = 5)
        @test_throws ArgumentError crop(h; start = 2.1, stop = 3)
        @test first(frequencies(crop(h; start = 6, stop = 6))) == 6
        shifted = FrequencySeries(a; frequencies = frequencies(h), epoch = t0 + 0.01, channels = names)
        @test_throws ArgumentError check_same_grid(h, shifted)
    end
    fine = TimeSeries(ones(4); epoch = t0, sampling_frequency = 1.0e9, channels = (:H1,))
    @test_throws ArgumentError crop(fine; start = t0, stop = t0)
    x = TimeSeries(ones(8, 2); epoch = t0, sampling_frequency = 10, channels = (:H1, :L1))
    @test_throws ArgumentError check_same_grid(x, TimeSeries(ones(8, 2); epoch = t0, sampling_frequency = 10, channels = (:L1, :H1)))
    @test_throws DimensionMismatch check_same_grid(x, TimeSeries(ones(7, 2); epoch = t0, sampling_frequency = 10, channels = (:H1, :L1)))
end

include("test_quality.jl")
