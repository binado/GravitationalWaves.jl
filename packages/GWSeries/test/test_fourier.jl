using FFTW

@testset "Physical Fourier convention" begin
    for T in (Float32, Float64), n in (2, 6, 16), nc in (1, 2)
        a = reshape(T.(sin.(1:(n * nc))), n, nc)
        nc == 1 && (a = vec(a))
        names = nc == 1 ? (:H1,) : (:H1, :L1)
        x = TimeSeries(a; epoch = 1126259460.125, sampling_frequency = 32, channels = names)
        h = to_frequencyseries(x)
        @test eltype(parent(h)) == Complex{T}
        @test parent(h) ≈ rfft(a, 1) ./ 32
        @test epoch(h) == epoch(x)
        @test channels(h) == names
        @test nsamples(h) == n ÷ 2 + 1
        y = to_timeseries(h)
        @test eltype(parent(y)) == T
        @test parent(y) ≈ a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        @test check_same_grid(x, y) === nothing
        shifted = TimeSeries(a; epoch = 0, sampling_frequency = 32, channels = names)
        @test parent(to_frequencyseries(shifted)) == parent(h)
    end
    x = TimeSeries(collect(1.0:16); epoch = 0, sampling_frequency = 16, channels = (:H1,))
    h = to_frequencyseries(x)
    power = abs2.(parent(h))
    @test sum(abs2, parent(x)) * sample_spacing(x) ≈ frequency_spacing(h) * (power[1] + power[end] + 2sum(power[2:(end - 1)]))
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
    precise = TimeSeries(Float32[1, 2, 3, 4]; epoch = big"1126259460.123456789", sampling_frequency = big"32", channels = (:H1,))
    @test epoch(to_timeseries(to_frequencyseries(precise))) == epoch(precise)
    @test parent(to_timeseries(to_frequencyseries(precise))) ≈ parent(precise)
    real_spectrum = FrequencySeries(ones(Float32, 5); frequencies = 0:1:4, epoch = 0, channels = (:H1,))
    @test parent(to_timeseries(real_spectrum)) ≈ Float32[8, 0, 0, 0, 0, 0, 0, 0]
    dirty = copy(h)
    parent(dirty)[1] += 7im
    parent(dirty)[end] -= 9im
    @test parent(to_timeseries(dirty)) == parent(to_timeseries(h))
    @test parent(to_frequencyseries(TimeSeries(view(parent(x), :); epoch = 0, sampling_frequency = 16, channels = (:H1,)))) ≈ parent(h)
end

@testset "Fourier conversion rejects unsupported storage with a clear error" begin
    x = TimeSeries(ones(Float64, 8, 3); epoch = 0, sampling_frequency = 8, channels = (:H1, :L1, :V1))
    h = to_frequencyseries(x)
    strided = TimeSeries(view(ones(16), 1:2:16); epoch = 0, sampling_frequency = 8, channels = (:H1,))
    @test nsamples(to_frequencyseries(strided)) == 5
    indexed = TimeSeries(view(parent(x), :, [1, 3]); epoch = 0, sampling_frequency = 8, channels = (:H1, :V1))
    for call in (
            () -> to_frequencyseries(indexed),
            () -> to_frequencyseries(TimeSeries(collect(1:8); epoch = 0, sampling_frequency = 8, channels = (:H1,))),
            () -> to_frequencyseries(TimeSeries(ones(ComplexF64, 8); epoch = 0, sampling_frequency = 8, channels = (:H1,))),
            () -> to_frequencyseries(h),
            () -> to_timeseries(x),
        )
        err = try
            call()
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("needs", err.msg) && occursin("got", err.msg)
    end
end

@testset "Analytic Fourier amplitudes" begin
    # A unit impulse at t = 0 has spectrum dt at every bin; a tone has amplitude N*dt/2 at its bin.
    n, fs = 64, 32.0
    impulse = TimeSeries([1.0; zeros(n - 1)]; epoch = 0, sampling_frequency = fs, channels = (:H1,))
    @test all(parent(to_frequencyseries(impulse)) .≈ 1 / fs)
    tone = TimeSeries(cos.(2π * 4 / n .* (0:(n - 1))); epoch = 0, sampling_frequency = fs, channels = (:H1,))
    spectrum = to_frequencyseries(tone)
    @test abs(parent(spectrum)[5]) ≈ n / fs / 2
    @test frequencies(spectrum)[5] ≈ 2.0
    @test sum(abs, parent(spectrum)[[1:4; 6:end]]) < 1.0e-12
end
