@testset "Series construction and metadata" begin
    a = reshape(collect(1.0:24.0), 8, 3)
    x = TimeSeries(a; epoch = 1126259460.125, sampling_frequency = 8, channels = ("H1", "L1", "V1"))
    @test x isa AbstractMatrix{Float64}
    @test x isa AbstractDimArray
    @test size(x) == (8, 3)
    @test parent(x) === a
    raw, meta = flatten(x)
    @test raw === a
    @test grid_metadata(TimeSeries(raw; meta...)) == meta
    @test step(x) == sampling_interval(x) == 0.125
    @test sampling_frequency(x) == 8
    @test duration(x) == 1
    # Accessors return raw Base ranges; the grid types live in the lookups.
    @test times(x) isa StepRangeLen
    @test parent(lookup(x, GPSTime)) isa TimeGrid
    @test parent(parent(lookup(x, GPSTime))) === times(x)
    @test times(x) == range(epoch(x); step = 0.125, length = 8)
    @test lookup(x, GPSTime) == times(x)
    @test channels(x) == (:H1, :L1, :V1)
    @test grid_metadata(TimeSeries(zeros(4, 3); meta...)) == meta
    v = ones(4)
    onechannel = TimeSeries(v; epoch = 0, sampling_interval = 1, channels = (:H1,))
    onechannel[2, 1] = 3
    @test v[2] == 3
    @test size(onechannel) == (4, 1)
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = 1, sampling_frequency = 1, channels = channels(x))
    @test_throws ArgumentError TimeSeries(a; epoch = Inf, sampling_interval = 1, channels = channels(x))
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = -1, channels = channels(x))
    @test_throws DimensionMismatch TimeSeries(a; epoch = 0, sampling_interval = 1, channels = (:H1,))
    @test_throws ArgumentError TimeSeries(a; epoch = 0, sampling_interval = 1, channels = (:H1, :H1, :V1))
    b = TimeSeries(zeros(Float32, 4); epoch = big"1126259460.123456789", sampling_interval = big"0.125", channels = (:H1,))
    @test epoch(b) isa BigFloat
    @test sampling_interval(b) isa BigFloat
    @test epoch(b) == big"1126259460.123456789"
    for n in (0, 1)
        s = TimeSeries(zeros(n); epoch = 2, sampling_interval = 0.25, channels = (:H1,))
        @test nsamples(s) == n
        @test epoch(s) == 2
    end
    # Grids are bits types, so series wrappers stay cheap to copy.
    @test isbitstype(typeof(parent(lookup(x, GPSTime))))
    h = FrequencySeries(complex.(a); epoch = 100, frequency_spacing = 0.5, channels = channels(x))
    @test isbitstype(typeof(parent(lookup(h, Freq))))
    @test parent(lookup(h, Freq)) isa SpectralGrid
    @test frequencies(h) === range(0.0; step = 0.5, length = 8)
end

@testset "Grid ranges" begin
    g = TimeGrid(1126259460.125, 0.25, 8)
    @test g == range(1126259460.125; step = 0.25, length = 8)
    @test g[3] == 1126259460.625
    @test_throws BoundsError g[9]
    @test g[2:8] isa TimeGrid && first(g[2:8]) == g[2] && step(g[2:8]) == 0.25
    @test g[1:3:8] isa TimeGrid && step(g[1:3:8]) == 0.75 && length(g[1:3:8]) == 3
    @test g[1:0] isa TimeGrid && isempty(g[1:0])
    @test !(g[8:-1:1] isa TimeGrid) && g[8:-1:1] == reverse(collect(g))
    @test !(reverse(g) isa TimeGrid) && reverse(g) == reverse(collect(g))
    @test view(g, 2:4) === g[2:4]
    # Values come from the wrapped Base range, so slices agree bitwise with it.
    @test parent(g) === range(1126259460.125; step = 0.25, length = 8)
    @test collect(g[2:8]) == collect(parent(g)[2:8])
    nondyadic = TimeGrid(-0.3, 0.1, 50)[2:49]
    @test collect(nondyadic) == collect(TimeGrid(first(nondyadic), step(nondyadic), length(nondyadic)))
    @test isfinite(step(TimeGrid(0.0, 1.0e300, 3)[1:typemax(Int):3]))
    f = SpectralGrid(0.5, 8, 100.0, 14)
    @test f == 0:0.5:3.5
    @test f[1:8] === f
    @test parent(f) === range(0.0; step = 0.5, length = 8)
    @test !(f[2:8] isa SpectralGrid) && f[2:8] == collect(f)[2:8]
    @test f[1:2:7] == collect(f)[1:2:7]
    c = GWSeries.ChannelNames((:H1, :L1, :V1))
    @test c[[3, 1]] == GWSeries.ChannelNames((:V1, :H1))
    @test c[2] === :L1
    @test Tuple(c) === (:H1, :L1, :V1)
end

using ForwardDiff
@testset "Grids with dual-number coordinates" begin
    # Waveform epochs can depend on differentiated parameters such as the total mass.
    function objective(m)
        x = TimeSeries(ones(4, 2); epoch = -m * 1.0e-3, sampling_interval = m * 1.0e-4, channels = (:plus, :cross))
        return epoch(x[2:3, :]) + sum(x[Chan = At(:plus)]) * sampling_interval(x)
    end
    @test ForwardDiff.derivative(objective, 60.0) ≈ -1.0e-3 + 1.0e-4 + 4 * 1.0e-4
end

@testset "Selection and slicing" begin
    a = reshape(collect(1.0:24.0), 8, 3)
    x = TimeSeries(a; epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1, :V1))
    @test x[2, 3] == a[2, 3]
    @test parent(x[Chan = At(:H1)]) == a[:, 1]
    @test x[Chan = At(:H1)] isa DimArray
    @test lookup(x[Chan = At(:H1)], GPSTime) == times(x)
    s = x[GPSTime = 2:8, Chan = At([:V1, :H1])]
    @test s isa TimeSeries
    @test channels(s) == (:V1, :H1)
    @test parent(s) == a[2:8, [3, 1]]
    @test epoch(s) == 100.25
    @test sampling_interval(s) == 0.25
    # Strided slices keep a uniform grid with a multiplied step.
    strided = x[2:2:8, :]
    @test strided isa TimeSeries
    @test epoch(strided) == 100.25 && sampling_interval(strided) == 0.5
    window = x[GPSTime = 100.5 .. 101.0]
    @test window isa TimeSeries
    @test epoch(window) == 100.5 && nsamples(window) == 3
    @test x[GPSTime = At(100.5), Chan = At(:L1)] == a[3, 2]
    @test_throws DimensionalData.Lookups.SelectorError x[Chan = At(:missing)]
    # Repeated channels and reordered samples keep labels but are no longer series.
    @test x[:, [1, 1]] isa DimArray && !(x[:, [1, 1]] isa AbstractSeries)
    @test !(x[[1, 3, 4], :] isa AbstractSeries)
    @test !(x[end:-1:1, :] isa AbstractSeries)
    @test x[1:0, :] isa TimeSeries
    v = view(x; GPSTime = 2:8, Chan = At([:L1, :H1]))
    @test v isa TimeSeries
    @test parent(v) isa SubArray
    v[1, 1] = 99
    @test x[2, 2] == 99
    h = FrequencySeries(complex.(a); epoch = 100, frequency_spacing = 0.5, channels = channels(x), time_samples = 14)
    # Partial bands are not complete spectra: they drop epoch and time length.
    for band in (h[3:2:7, :], h[2:4, :], h[Freq = 1 .. 2], view(h, 2:4, :))
        @test band isa DimArray && !(band isa AbstractSeries)
    end
    @test parent(h[Freq = 1 .. 2]) == parent(h)[3:5, :]
    @test lookup(h[Freq = 1 .. 2], Freq) == 1.0:0.5:2.0
    full = h[Chan = At([:L1])]
    @test full isa FrequencySeries
    @test time_samples(full) == 14 && epoch(full) == 100
    @test channels(view(h; Chan = At([:L1, :H1]))) == (:L1, :H1)
    @test frequencies(h) == 0:0.5:3.5
    @test frequency_bins(h; fmin = 0.6, fmax = 1.6) == 3:4
    @test frequency_bins(h; fmin = -1, fmax = 1) == 1:3
    @test frequency_bins(h) == axes(parent(h), 1)
    @test isempty(frequency_bins(h; fmin = 4, fmax = 5))
    @test isempty(frequency_bins(h; fmin = -2, fmax = -1))
    @test_throws ArgumentError frequency_bins(h; fmin = 2, fmax = 1)
    @test_throws ArgumentError frequency_bins(h; fmax = Inf)
    for kw in ((; frequency_start = 0), (; frequencies = 0:0.5:3.5), (; spectrum = :onesided))
        @test_throws MethodError FrequencySeries(parent(h); grid_metadata(h)..., kw...)
    end
end

@testset "Metadata-preserving operations" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    b = ones(8, 2)
    for y in (
            2 .* x, x .+ b, b .+ x, x + b, b + x, x - b, x + x, 2x, x / 2, -x, map(sin, x), map(+, x, x),
            x .* x, copy(x), similar(x), similar(x, Float32), cumsum(x; dims = 1), cumsum(x; dims = GPSTime),
        )
        @test y isa TimeSeries
        @test grid_metadata(y) == grid_metadata(x)
    end
    @test eltype(similar(x, Float32)) === Float32
    @test parent(@. 2 * x + b) == 2 .* parent(x) .+ b
    # Products, reorderings, reductions and channel accumulation have no series grid.
    mix = [0.0 1.0; 1.0 0.0]
    for y in (
            x * mix, x * x', x * ones(2, 3), permutedims(x), reverse(x; dims = 1), circshift(x, 1),
            diff(x; dims = 1), sum(x; dims = 1), cumsum(x; dims = 2), maximum(abs, x; dims = 1),
        )
        @test y isa AbstractDimArray && !(y isa AbstractSeries)
    end
    @test parent(x * mix) == parent(x) * mix
    @test parent(diff(x; dims = 1)) == diff(parent(x); dims = 1)
    @test parent(reverse(x; dims = 1)) == reverse(parent(x); dims = 1)
    @test parent(circshift(x, 1)) == circshift(parent(x), 1)
    @test sum(x) == sum(parent(x))
    out = copy(x)
    @test out !== x
    out .= 0
    @test all(iszero, parent(out))
    out .= 2 .* x .+ b
    @test parent(out) == 2 .* parent(x) .+ b
    raw = similar(parent(x))
    raw .= 2 .* x
    @test raw == 2 .* parent(x)
    shifted = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = 101)
    @test_throws ArgumentError x .+ shifted
    @test_throws ArgumentError out .= shifted
    @test_throws ArgumentError x + shifted
    @test_throws ArgumentError map(+, x, shifted)
    @test_throws ArgumentError copyto!(out, shifted)
    swapped = TimeSeries(copy(parent(x)); grid_metadata(x)..., channels = (:L1, :H1))
    @test_throws DimensionMismatch x .+ swapped
    @test_throws ArgumentError check_same_grid(x, swapped)
    close = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = nextfloat(epoch(x)))
    @test check_same_grid(x, close) === nothing
    @test x .+ close isa TimeSeries
    @test check_same_grid(x, x) === nothing
end

@testset "Epoch shifts beyond DimensionalData's relative tolerance throw" begin
    # DimensionalData alone compares lookups with a relative ≈, which accepts GPS epochs
    # seconds apart; the GPSTime comparison hook must reject them. This guards DD upgrades.
    t0 = 1126259460.125
    x = TimeSeries(randn(MersenneTwister(1), 64, 1); epoch = t0, sampling_frequency = 64, channels = (:H1,))
    one_second = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = t0 + 1)
    one_sample = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = t0 + 1 / 64)
    @test isapprox(t0, t0 + 1)
    for y in (one_second, one_sample)
        @test_throws "epochs differ" x .+ y
        @test_throws "epochs differ" x + y
    end
    @test x .+ TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = nextfloat(t0)) isa TimeSeries
    h = rfft(x)
    later = FrequencySeries(copy(parent(h)); grid_metadata(h)..., epoch = t0 + 1)
    @test_throws "epochs differ" h .+ later
    # The predicate form of the comparison reports a mismatch instead of throwing.
    @test !DimensionalData.Dimensions.comparedims(Bool, x, one_second; val = true)
    @test DimensionalData.Dimensions.comparedims(Bool, x, x; val = true)
end

@testset "Crop and pad" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    y = crop(x; start = 100.5, stop = 101.5)
    @test y isa TimeSeries
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

struct ZeroBasedMatrix <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end
Base.size(x::ZeroBasedMatrix) = size(x.data)
Base.axes(x::ZeroBasedMatrix) = (0:(size(x, 1) - 1), axes(x.data, 2))
Base.getindex(x::ZeroBasedMatrix, i::Int, j::Int) = x.data[i + 1, j]

@testset "Faithful dimensional reconstruction" begin
    DD = DimensionalData
    L = DD.Lookups
    x = TimeSeries(ones(8, 2); epoch = 1126259460.125, sampling_interval = 0.25, channels = (:H1, :L1))
    for s in (x, rfft(x))
        d = dims(s)
        l = lookup(s, 1)
        # Canonical lookups, including their metadata, are retained verbatim.
        tagged = rebuild(l; metadata = Dict(:source => "measurement"))
        td = (rebuild(d[1], tagged), d[2])
        canonical = rebuild(s; dims = td)
        @test canonical isa (s isa TimeSeries ? TimeSeries : FrequencySeries)
        @test parent(canonical) === parent(s)
        @test dims(canonical) === td
        @test DD.metadata(lookup(canonical, 1)) === DD.metadata(tagged)
        for modified in (
                rebuild(l; sampling = L.Intervals(L.Start())),
                rebuild(l; order = L.ReverseOrdered()),
                rebuild(l; span = L.Irregular((nothing, nothing))),
                rebuild(l; span = L.Regular(nextfloat(step(parent(l))))),
                rebuild(l; sampling = L.Intervals(L.Center()), metadata = Dict(:unit => :seconds)),
            )
            supplied = (rebuild(d[1], modified), d[2])
            y = rebuild(s; dims = supplied)
            @test y isa DimArray
            @test dims(y) === supplied
            @test lookup(y, 1) === modified
            @test parent(y) === parent(s)
        end
        @test_throws DimensionMismatch rebuild(s; data = zeros(1, 2))
        @test_throws DimensionMismatch rebuild(s; data = ZeroBasedMatrix(ones(size(s))))
        @test_throws DimensionMismatch rebuild(s; data = zeros(size(s, 1), 1))
        @test dims(copy(s)) === d
        @test parent(copy(s)) !== parent(s)
        @test parent(view(s, :, :)) isa SubArray
        @test !(s[:, [1, 1]] isa AbstractSeries)
        # Canonical series deliberately carry no DD array name, refdims, metadata.
        renamed = rebuild(s; name = :ignored, refdims = (DD.Z(1),), metadata = Dict(:a => 1))
        @test DD.name(renamed) isa DD.NoName
        @test DD.refdims(renamed) == ()
        @test DD.metadata(renamed) isa DD.NoMetadata
    end
    @test x[:, 1:0] isa TimeSeries
    @test x[1:0, :] isa TimeSeries
    @test x[1:2:end, :] isa TimeSeries
end

@testset "Mixed labelled operations validate before mutation" begin
    DD = DimensionalData
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 1126259460.125, sampling_interval = 0.25, channels = (:H1, :L1))
    labelled = DimArray(copy(parent(x)), dims(x))
    @test map(+, x, labelled) isa TimeSeries
    @test map(+, labelled, x) isa DimArray
    @test parent(map(+, labelled, x)) == 2 .* parent(x)
    @test parent(map(+, x, labelled, x)) == 3 .* parent(x)
    for other in (
                TimeSeries(fill(-1.0, size(x)); grid_metadata(x)..., epoch = epoch(x) + 1),
                TimeSeries(fill(-1.0, size(x)); grid_metadata(x)..., sampling_interval = 0.5),
                TimeSeries(fill(-1.0, size(x)); grid_metadata(x)..., channels = (:L1, :H1)),
            ), source in (
                other, DimArray(parent(other), dims(other)),
                DimArray(parent(other), (GPSTime(times(other)), dims(other, Chan))),
                DimArray(parent(other), (GPSTime(collect(times(other))), dims(other, Chan))),
            )
        @test_throws Union{ArgumentError, DimensionMismatch} map(+, x, source)
        @test_throws Union{ArgumentError, DimensionMismatch} map(+, source, x)
        @test_throws Union{ArgumentError, DimensionMismatch} map(+, x, labelled, source)
        out = copy(x)
        saved = copy(parent(out))
        @test_throws Union{ArgumentError, DimensionMismatch} copyto!(out, source)
        @test parent(out) == saved
        @test_throws Union{ArgumentError, DimensionMismatch} out[:, :] = source
        @test parent(out) == saved
    end
    for coordinates in (times(x), collect(times(x)))
        plain = DimArray(copy(parent(x)), (GPSTime(coordinates), dims(x, Chan)))
        @test parent(map(+, plain, x)) == 2 .* parent(x)
        @test parent(map(+, x, plain)) == 2 .* parent(x)
        @test parent(copyto!(copy(x), plain)) == parent(x)
    end
    out = copy(x)
    @test copyto!(out, labelled) === out
    out[GPSTime = 2:4, Chan = At([:L1])] = DimArray(fill(42.0, 3, 1), dims(x[2:4, [2]]))
    @test all(==(42), parent(out)[2:4, 2])
    saved = copy(parent(out))
    @test_throws ArgumentError out[2:4, :] = DimArray(parent(x)[1:3, :], dims(x[1:3, :]))
    @test parent(out) == saved
    out[:, :] = zeros(size(x))
    @test all(iszero, parent(out))
    @test copyto!(out, ones(size(x))) === out
    @test all(isone, parent(out))
    # Empty indices must reach the terminal validation without dispatch ambiguity.
    out[Int[], :] = DimArray(zeros(0, 2), dims(x[1:0, :]))
    h = rfft(x)
    odd = FrequencySeries(fill(0.0im, size(h)); grid_metadata(h)..., time_samples = 9)
    for source in (odd, DimArray(parent(odd), dims(odd)))
        @test_throws "original time lengths differ" map(+, h, source)
        @test_throws "original time lengths differ" map(+, source, h)
        out = copy(h)
        saved = copy(parent(out))
        @test_throws "original time lengths differ" copyto!(out, source)
        @test_throws "original time lengths differ" out[:, :] = source
        @test parent(out) == saved
    end
    @test DD.Dimensions.comparedims(Bool, x, x; val = true)
    @test parent(x .+ x[:, [1]]) == parent(x) .+ parent(x)[:, [1]]
    singleton = TimeSeries(ones(1, 2); grid_metadata(x)...)
    @test parent(x .+ singleton) == parent(x) .+ 1
end
