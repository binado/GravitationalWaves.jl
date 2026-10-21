@testset "Series construction and metadata" begin
    a = reshape(collect(1.0:24.0), 8, 3)
    x = TimeSeries(a; epoch = 1126259460.125, sampling_frequency = 8, channels = ("H1", "L1", "V1"))
    @test x isa AbstractMatrix{Float64}
    @test size(x) == (8, 3)
    @test parent(x) === a
    raw, meta = flatten(x)
    @test raw === a
    @test metadata(TimeSeries(raw; meta...)) == meta
    @test step(x) == sampling_interval(x) == 0.125
    @test sampling_frequency(x) == 8
    @test duration(x) == 1
    @test times(x) == range(epoch(x); step = 0.125, length = 8)
    @test channels(x) == (:H1, :L1, :V1)
    @test metadata(TimeSeries(zeros(4, 3); meta...)) == meta
    v = ones(4)
    onechannel = TimeSeries(v; epoch = 0, sampling_interval = 1, channels = (:H1,))
    onechannel[2, 1] = 3
    @test v[2] == 3
    @test size(onechannel) == (4, 1)
    @test setindex!(onechannel, 4, 1, 1) === onechannel
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = 1, sampling_frequency = 1, channels = channels(x))
    @test_throws ArgumentError TimeSeries(a; epoch = Inf, sampling_interval = 1, channels = channels(x))
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = -1, channels = channels(x))
    @test_throws DimensionMismatch TimeSeries(a; epoch = 0, sampling_interval = 1, channels = (:H1,))
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = 1, channels = (:H1, :H1, :V1))
    b = TimeSeries(zeros(Float32, 4); epoch = big"1126259460.123456789", sampling_interval = big"0.125", channels = (:H1,))
    @test epoch(b) isa BigFloat
    @test sampling_interval(b) isa BigFloat
    for n in (0, 1)
        s = TimeSeries(zeros(n); epoch = 2, sampling_interval = 0.25, channels = (:H1,))
        @test nsamples(s) == n
        @test first(times(s)) == 2
    end
end

@testset "Channel and sample indexing" begin
    a = reshape(collect(1.0:24.0), 8, 3)
    x = TimeSeries(a; epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1, :V1))
    @test x[2, 3] == a[2, 3]
    @test x[:, :H1] == a[:, 1]
    @test x[:, :H1] isa Vector
    @test x[2, :] isa Vector
    s = x[2:8, [:V1, :H1]]
    @test s isa TimeSeries
    @test channels(s) == (:V1, :H1)
    @test parent(s) == a[2:8, [3, 1]]
    @test epoch(s) == 100.25
    @test sampling_interval(s) == 0.25
    @test size(x[:, (:L1,)]) == (8, 1)
    @test_throws ArgumentError x[:, :missing]
    @test x[:, [1, 1]] isa Matrix
    @test x[[1, 3, 4], :] isa Matrix
    @test x[end:-1:1, :] isa Matrix
    @test x[1:0, :] isa TimeSeries
    v = view(x, 2:8, [:L1, :H1])
    @test v isa TimeSeries
    @test parent(v) isa SubArray
    v[1, 1] = 99
    @test x[2, :L1] == 99
    @test x[2:4, :].data !== a
    h = FrequencySeries(complex.(a); epoch = 100, frequency_spacing = 0.5, channels = channels(x))
    @test h[3:2:7, :] isa Matrix
    @test x[2:2:8, :] isa Matrix
    @test h[axes(h, 1), [:L1]] isa FrequencySeries
    @test view(h, 2:4, :) isa SubArray
    @test channels(view(h, :, [:L1, :H1])) == (:L1, :H1)
    @test h[:, :H1] isa Vector
    @test h[:, Symbol[]] isa Matrix
    @test h[:, [:H1, :H1]] isa Matrix
    @test frequencies(h) == 0:0.5:3.5
    @test frequency_bins(h; fmin = 0.6, fmax = 1.6) == 3:4
    @test frequency_bins(h; fmin = -1, fmax = 1) == 1:3
    @test frequency_bins(h) == axes(h, 1)
    @test isempty(frequency_bins(h; fmin = 4, fmax = 5))
    @test isempty(frequency_bins(h; fmin = -2, fmax = -1))
    @test_throws ArgumentError frequency_bins(h; fmin = 2, fmax = 1)
    @test_throws ArgumentError frequency_bins(h; fmax = Inf)
    for kw in ((; frequency_start = 0), (; frequencies = 0:0.5:3.5), (; spectrum = :onesided))
        @test_throws MethodError FrequencySeries(parent(h); metadata(h)..., kw...)
    end
    huge = TimeSeries(ones(1); epoch = 0, sampling_interval = 1.0e300, channels = (:H1,))
    @test huge[1:typemax(Int):1, :] isa Matrix
end

@testset "Metadata-preserving operations" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    b = ones(8, 2)
    for y in (2 .* x, x .+ b, b .+ x, x + b, b + x, x - b, 2x, x / 2, -x, map(sin, x), x .* x)
        @test y isa TimeSeries
        @test metadata(y) == metadata(x)
    end
    @test parent(@. 2 * x + b) == 2 .* parent(x) .+ b
    # Matrix products carry no grid: channel-mixing must not keep the channel labels.
    mix = [0.0 1.0; 1.0 0.0]
    @test x * mix isa Matrix
    @test x * mix == parent(x) * mix
    @test x * x' isa Matrix
    @test x * ones(2, 3) isa Matrix
    # Reordering or accumulating Base routines must not keep the source grid.
    for y in (reverse(x; dims = 1), circshift(x, 1), cumsum(x; dims = 1), diff(x; dims = 1), similar(x), similar(x, Float32))
        @test y isa Matrix
    end
    @test diff(x; dims = 1) == diff(parent(x); dims = 1)
    @test reverse(x; dims = 1) == reverse(parent(x); dims = 1)
    singleton = TimeSeries(ones(1, 2); epoch = 0, sampling_interval = 1, channels = (:H1, :L1))
    @test singleton .+ b isa Matrix
    @test reshape(x, :) isa Vector
    @test transpose(x) isa Transpose
    @test sum(x; dims = 1) isa Matrix
    @test maximum(abs, x; dims = 1) isa Matrix
    @test minimum(abs, x; dims = 2) isa Matrix
    @test prod(abs, x; dims = 1) isa Matrix
    @test similar(x, (4, 2)) isa Matrix
    @test copy(x) !== x
    out = copy(x)
    out .= 0
    @test all(iszero, parent(out))
    out .= 2 .* x .+ b
    @test parent(out) == 2 .* parent(x) .+ b
    raw = similar(parent(x))
    raw .= 2 .* x
    @test raw == 2 .* parent(x)
    shifted = TimeSeries(copy(parent(x)); metadata(x)..., epoch = 101)
    @test_throws ArgumentError x .+ shifted
    @test_throws ArgumentError out .= shifted
    @test_throws ArgumentError x + shifted
    @test_throws ArgumentError map(+, x, shifted)
    @test_throws ArgumentError map!(identity, out, shifted)
    @test_throws ArgumentError copyto!(out, shifted)
    close = TimeSeries(copy(parent(x)); metadata(x)..., epoch = nextfloat(epoch(x)))
    @test check_same_grid(x, close) === nothing
    @test check_same_grid(x, x) === nothing
end

@testset "Crop and pad" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    y = crop(x; start = 100.5, stop = 101.5)
    @test parent(y) == parent(x)[3:6, :]
    @test epoch(y) == 100.5
    @test_throws ArgumentError crop(x; start = 100.1, stop = 101)
    z = pad(x; before = 2, after = 3)
    @test epoch(z) == 99.5
    @test size(z) == (13, 2)
    @test parent(z)[3:10, :] == parent(x)
    @test all(iszero, parent(z)[1:2, :])
    @test_throws ArgumentError pad(x; before = -1)
    f = FrequencySeries(ones(4); epoch = 0, frequency_spacing = 1, channels = (:H1,))
    @test_throws MethodError pad(f; before = 1)
    @test_throws MethodError crop(f; start = 0, stop = 1)
end

include("prototypes.jl")
@testset "Grid compatibility tolerance" begin
    t0 = 1126259460.1
    x = time_prototype(; nsamples = 8, epoch = t0, sampling_frequency = 16, channels = (:H1,))
    mismatch(; epoch = t0, fs = 16, n = 8) =
        time_prototype(; nsamples = n, epoch, sampling_frequency = fs, channels = (:H1,))
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

    h = frequency_prototype(; frequencies = 0:0.5:4, epoch = t0, channels = (:H1,))
    other(; f = 0:0.5:4, epoch = t0) = frequency_prototype(; frequencies = f, epoch, channels = (:H1,))
    @test check_same_grid(h, other(; f = range(0.0; step = 0.5, length = 9))) === nothing
    @test check_same_grid(h, other(; f = LinRange(0, 4, 9))) === nothing
    @test check_same_grid(h, other(; epoch = nextfloat(t0))) === nothing
    @test_throws ArgumentError check_same_grid(h, other(; f = 0:0.25:2))
    @test_throws ArgumentError check_same_frequency_grid(frequencies(h), 0.5:0.5:4.5)
    @test_throws ArgumentError check_same_grid(h, other(; f = 0:0.5001:4.0008))
    @test_throws ArgumentError check_same_grid(h, other(; epoch = t0 + 0.01))
    @test check_same_frequency_grid(0:0.5:4, LinRange(0, 4, 9)) === nothing
    @test_throws DimensionMismatch check_same_frequency_grid(0:0.5:4, 0:0.5:3.5)
    @test_throws ArgumentError check_same_frequency_grid(0:0.5:4, 0:0.25:2)
end

@testset "Identical coordinates at high rates" begin
    gps = 1.2e9
    fast = time_prototype(; nsamples = 4, epoch = gps, sampling_frequency = 1.0e6, channels = (:H1,))
    @test check_same_grid(fast, fast) === nothing
    @test nsamples(crop(fast; start = gps, stop = gps)) == 0
    err = try
        crop(fast; start = gps, stop = gps + 2.0e-6)
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
    d = time_prototype(; epoch = 0, nsamples = 4, sampling_frequency = 4, channels = (:H1,))
    @test_throws ArgumentError check_data(d, ZeroBasedVector(ones(4)))
    @test check_data(d, view(ones(8), 1:2:8)) === nothing
    e = frequency_prototype(; frequencies = range(0.0; step = 0.5, length = 0), epoch = 0, channels = (:H1,))
    @test frequency_spacing(e) == 0.5
    @test isempty(frequency_bins(e))
    @test_throws MethodError pad(e; after = 2)
end
