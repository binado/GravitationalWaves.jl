using Test
using GWTimeFrequencyDomains

@testset "Detached domains and data checks" begin
    td = TimeDomain(; epoch = 1126259460.125, nsamples = 8, sampling_frequency = 16, channels = ("H1",))
    @test !(td isa AbstractArray)
    @test channels(td) == (:H1,)
    @test nsamples(td) == 8
    @test sample_spacing(td) == 1 / 16
    @test sampling_frequency(td) == 16
    @test duration(td) == 0.5
    @test times(td) == range(epoch(td); step = 1 / 16, length = 8)
    @test check_data(td, ones(8)) === nothing
    @test check_data(td, ones(8, 1)) === nothing
    @test_throws DimensionMismatch check_data(td, ones(7))
    @test_throws DimensionMismatch check_data(td, ones(8, 2))
    @test_throws ArgumentError check_data(td, ["a"])
    @test_throws ArgumentError check_data(td, ones(8, 1, 1))
    @test_throws ArgumentError TimeDomain(; epoch = 0, nsamples = -1, sampling_frequency = 16, channels = (:H1,))
    for fs in (0, -1, Inf, NaN)
        @test_throws ArgumentError TimeDomain(; epoch = 0, nsamples = 8, sampling_frequency = fs, channels = (:H1,))
    end
    for names in ((), (:H1, "H1"), (1,))
        @test_throws ArgumentError TimeDomain(; epoch = 0, nsamples = 8, sampling_frequency = 16, channels = names)
    end
    @test_throws ArgumentError TimeDomain(; epoch = Inf, nsamples = 8, sampling_frequency = 16, channels = (:H1,))
    big = TimeDomain(; epoch = big"1126259460.123456789", nsamples = 8, sampling_frequency = big"16", channels = (:H1,))
    @test epoch(big) isa BigFloat
    @test sample_spacing(big) isa BigFloat
    fd = frequency_domain(td)
    @test nsamples(fd) == 5
    @test frequency_spacing(fd) == 2
    @test frequencies(fd) == 0:2:8
    @test check_same_grid(td, time_domain(fd)) === nothing
    @test_throws MethodError duration(fd)
    @test_throws ArgumentError check_same_grid(td, fd)
    for f in (-1:1:2, 4:-1:1, range(0.0; step = Inf, length = 4))
        @test_throws ArgumentError FrequencyDomain(; frequencies = f, epoch = 0, channels = (:H1,))
    end
    for n in (0, 1, 3)
        @test_throws ArgumentError check_fourier_domain(TimeDomain(; epoch = 0, nsamples = n, sampling_frequency = 16, channels = (:H1,)))
    end
    for f in (range(0.0; step = 1, length = 0), range(0.0; step = 1, length = 1), 1:1:4)
        @test_throws ArgumentError check_fourier_domain(FrequencyDomain(; frequencies = f, epoch = 0, channels = (:H1,)))
    end
end

@testset "Grid operations" begin
    t0 = 1126259460.1
    for a in (collect(1.0:8), reshape(collect(1.0:16), 8, 2))
        names = ndims(a) == 1 ? (:H1,) : (:H1, :L1)
        td = TimeDomain(; epoch = t0, nsamples = 8, sampling_frequency = 10, channels = names)
        d, data = crop(td, a; start = t0 + 0.2, stop = t0 + 0.6)
        @test nsamples(d) == 4
        @test epoch(d) == t0 + 0.2
        @test data == (ndims(a) == 1 ? a[3:6] : a[3:6, :])
        @test data !== a
        @test_throws ArgumentError crop(td, a; start = t0 + 0.21, stop = t0 + 0.6)
        @test_throws ArgumentError crop(td, a; start = t0 - 0.1, stop = t0)
        @test_throws ArgumentError crop(td, a; start = t0 + 0.6, stop = t0 + 0.2)
        empty, data = crop(td, a; start = t0 + 0.8, stop = t0 + 0.8)
        @test nsamples(empty) == 0 && isempty(data)
        @test epoch(empty) == t0 + 0.8
        padded, values = pad(td, a; before = 2, after = 3)
        @test nsamples(padded) == 13
        @test epoch(padded) == t0 - 0.2
        restored, data = crop(padded, values; start = t0, stop = t0 + 0.8)
        @test check_same_grid(restored, td) === nothing
        @test data == a
        @test_throws ArgumentError pad(td, a; before = -1)
        fd = FrequencyDomain(; frequencies = 2:0.5:5.5, epoch = t0, channels = names)
        padded, values = pad(fd, a; before = 4, after = 2)
        @test first(frequencies(padded)) == 0 && epoch(padded) == t0
        @test last(crop(padded, values; start = 2, stop = 6)) == a
        @test_throws ArgumentError pad(fd, a; before = 5)
        empty, values = crop(fd, a; start = 6, stop = 6)
        @test first(frequencies(empty)) == 6 && isempty(values)
    end
end

@testset "Grid compatibility tolerance" begin
    t0 = 1126259460.1
    x = TimeDomain(; nsamples = 8, epoch = t0, sampling_frequency = 16, channels = (:H1,))
    mismatch(; epoch = t0, fs = 16, n = 8) =
        TimeDomain(; nsamples = n, epoch, sampling_frequency = fs, channels = (:H1,))
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

    h = FrequencyDomain(; frequencies = 0:0.5:4, epoch = t0, channels = (:H1,))
    other(; f = 0:0.5:4, epoch = t0) = FrequencyDomain(; frequencies = f, epoch, channels = (:H1,))
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
    fast = TimeDomain(; nsamples = 4, epoch = gps, sampling_frequency = 1.0e6, channels = (:H1,))
    @test check_same_grid(fast, fast) === nothing
    @test nsamples(first(crop(fast, ones(4); start = gps, stop = gps))) == 0
    err = try
        crop(fast, ones(4); start = gps, stop = gps + 2.0e-6)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("crop boundary", err.msg) && occursin("spacing", err.msg)
end


struct ZeroBasedVector <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(x::ZeroBasedVector) = size(x.data)
Base.axes(x::ZeroBasedVector) = (0:(length(x.data) - 1),)
Base.getindex(x::ZeroBasedVector, i::Int) = x.data[i + 1]
@testset "Axes and views" begin
    d = TimeDomain(; epoch = 0, nsamples = 4, sampling_frequency = 4, channels = (:H1,))
    @test_throws ArgumentError check_data(d, ZeroBasedVector(ones(4)))
    @test check_data(d, view(ones(8), 1:2:8)) === nothing
    e = FrequencyDomain(; frequencies = range(0.0; step = 0.5, length = 0), epoch = 0, channels = (:H1,))
    @test frequency_spacing(e) == 0.5
    p, values = pad(e, Float64[]; after = 2)
    @test nsamples(p) == 2 && values == [0, 0]
end
include("test_windows.jl")
include("test_fourier.jl")
include("test_quality.jl")
