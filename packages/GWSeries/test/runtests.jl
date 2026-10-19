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
    @test nsamples(crop(fine; start = t0, stop = t0)) == 0
    @test_throws ArgumentError crop(fine; start = t0, stop = t0 + 1.0e-6)
    x = TimeSeries(ones(8, 2); epoch = t0, sampling_frequency = 10, channels = (:H1, :L1))
    @test_throws ArgumentError check_same_grid(x, TimeSeries(ones(8, 2); epoch = t0, sampling_frequency = 10, channels = (:L1, :H1)))
    @test_throws DimensionMismatch check_same_grid(x, TimeSeries(ones(7, 2); epoch = t0, sampling_frequency = 10, channels = (:H1, :L1)))
end

@testset "Grid compatibility tolerance" begin
    t0 = 1126259460.1
    x = TimeSeries(ones(8); epoch = t0, sampling_frequency = 16, channels = (:H1,))
    mismatch(; epoch = t0, fs = 16, n = 8) =
        TimeSeries(ones(n); epoch, sampling_frequency = fs, channels = (:H1,))
    @test check_same_grid(x, mismatch(; epoch = nextfloat(t0))) === nothing
    @test check_same_grid(x, mismatch(; epoch = prevfloat(t0, 2))) === nothing
    @test_throws ArgumentError check_same_grid(x, mismatch(; epoch = t0 + 1.0e-3))
    # Documented resolution: offsets of a few ulp of a GPS epoch are the same epoch.
    @test check_same_grid(x, mismatch(; epoch = t0 + 1.0e-6)) === nothing
    @test_throws ArgumentError check_same_grid(x, mismatch(; epoch = t0 + 1.0e-5))
    @test_throws ArgumentError check_same_grid(x, mismatch(; fs = 32))
    @test check_same_grid(x, mismatch(; fs = nextfloat(16.0))) === nothing
    @test_throws DimensionMismatch check_same_grid(x, mismatch(; n = 7))
    err = try
        check_same_grid(x, mismatch(; epoch = t0 + 0.5))
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("epochs differ", err.msg) && occursin("same start time", err.msg)

    h = FrequencySeries(ones(9); frequencies = 0:0.5:4, epoch = t0, channels = (:H1,))
    other(; f = 0:0.5:4, epoch = t0) = FrequencySeries(ones(length(f)); frequencies = f, epoch, channels = (:H1,))
    @test check_same_grid(h, other(; f = range(0.0; step = 0.5, length = 9))) === nothing
    @test check_same_grid(h, other(; f = LinRange(0, 4, 9))) === nothing
    @test check_same_grid(h, other(; epoch = nextfloat(t0))) === nothing
    @test_throws ArgumentError check_same_grid(h, other(; f = 0:0.25:2))
    @test_throws ArgumentError check_same_grid(h, other(; f = 0.5:0.5:4.5))
    @test_throws ArgumentError check_same_grid(h, other(; f = 0:0.5001:4.0008))
    @test_throws ArgumentError check_same_grid(h, other(; epoch = t0 + 0.01))
    @test check_same_frequency_grid(0:0.5:4, LinRange(0, 4, 9)) === nothing
    @test_throws DimensionMismatch check_same_frequency_grid(0:0.5:4, 0:0.5:3.5)
    @test_throws ArgumentError check_same_frequency_grid(0:0.5:4, 0:0.25:2)
end

@testset "Identical coordinates at high rates" begin
    gps = 1.2e9
    fast = TimeSeries(ones(4); epoch = gps, sampling_frequency = 1.0e6, channels = (:H1,))
    @test check_same_grid(fast, copy(fast)) === nothing
    @test nsamples(crop(fast; start = gps, stop = gps)) == 0
    err = try
        crop(fast; start = gps, stop = gps + 2.0e-6)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("crop boundary", err.msg) && occursin("spacing", err.msg)
end

@testset "Internal constructor and parameter constraints" begin
    @test_throws TypeError TimeSeries{String, Float64, Float64, Tuple{Symbol}}
    @test_throws TypeError FrequencySeries{Vector{Float64}, Float64, Vector{Float64}, Tuple{Symbol}}
    @test_throws MethodError TimeSeries(ones(4), 0.0, 0.25, (:H1,), Val(:validated))
end

include("test_fourier.jl")
include("test_quality.jl")

# A backend-like array with unsupported axes; no extra runtime dependency needed.
struct ZeroBasedVector <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(x::ZeroBasedVector) = size(x.data)
Base.axes(x::ZeroBasedVector) = (0:(length(x.data) - 1),)
Base.getindex(x::ZeroBasedVector, i::Int) = x.data[i + 1]
@testset "Unsupported axes and storage views" begin
    a = ZeroBasedVector(ones(4))
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_frequency = 4, channels = (:H1,))
    @test_throws ArgumentError FrequencySeries(a; frequencies = 0:1:3, epoch = 0, channels = (:H1,))
    v = view(ones(8), 1:2:8)
    x = TimeSeries(v; epoch = 0, sampling_frequency = 4, channels = (:H1,))
    @test parent(x) === v
    @test_throws ArgumentError with_data(x, a)
    f = range(0.0; step = 0.5, length = 0)
    empty = FrequencySeries(Float64[]; frequencies = f, epoch = 0, channels = (:H1,))
    @test nsamples(empty) == 0
    @test first(frequencies(empty)) == 0
    @test frequency_spacing(empty) == 0.5
    @test nsamples(pad(empty; after = 2)) == 2
    @test check_same_grid(empty, copy(empty)) === nothing
end
