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
        @test grid_metadata(irfft(h, n)) == grid_metadata(x)
        @test parent(to_timedomain(to_frequencydomain(x); nsamples = n)) ≈ a
        @test_throws ArgumentError rfft(x, 2)
        @test_throws ArgumentError ifft(h)
        @test_throws ArgumentError fft(x)
        @test_throws ArgumentError fft!(x)
        @test_throws MethodError to_timedomain(h[2:end, :]; nsamples = n)
        @test h[2:end, :] isa DimArray
        @test_throws ArgumentError rfft(h)
        @test_throws ArgumentError irfft(x, n)
        @test_throws ArgumentError plan_rfft(x)
        @test time_samples(h) == n
        @test parent(to_timedomain(h)) ≈ a
        @test parent(irfft(h)) ≈ a
        @test parent(transform(plan_to_timedomain(h), h)) ≈ a
        # N and N+1 (or N-1) can share a bin count; the recorded length disambiguates.
        @test_throws ArgumentError irfft(h, n + 1)
        @test_throws ArgumentError irfft(h, n - 1)
        @test_throws ArgumentError to_timedomain(h; nsamples = n + 1)
        @test_throws ArgumentError irfft(h, 2n)
        let
            p = plan_to_frequencydomain(x)
            actual = transform(p, x)
            expected = h
            @test parent(actual) ≈ parent(expected)
            @test grid_metadata(p) == grid_metadata(expected)
            raw = transform(p, parent(x))
            @test raw ≈ parent(expected)
            target = copy(expected)
            transform!(target, p, x)
            @test parent(target) ≈ parent(expected)
            @test transform!(parent(target), p, x) === parent(target)
            @test transform!(target, p, parent(x)) === target
            wrong_epoch = FrequencySeries(copy(parent(target)); grid_metadata(target)..., epoch = epoch(target) + 1)
            @test_throws ArgumentError transform!(wrong_epoch, p, parent(x))
            shifted_input = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = epoch(x) + 1)
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

@testset "Recorded time length" begin
    x = TimeSeries(randn(MersenneTwister(5), 8); epoch = 0, sampling_interval = 0.25, channels = (:H1,))
    h = rfft(x)
    # Spectra that did not come from a series FFT leave the length unrecorded and need it explicitly.
    bare = FrequencySeries(copy(parent(h)); epoch = 0, frequency_spacing = frequency_spacing(h), channels = (:H1,))
    @test time_samples(bare) === nothing
    @test_throws ArgumentError irfft(bare)
    @test_throws ArgumentError to_timedomain(bare)
    @test irfft(bare, 8) isa TimeSeries
    @test irfft(bare, 9) isa TimeSeries
    # Recorded and unrecorded spectra share a grid; two different recorded lengths do not.
    @test check_same_grid(h, bare) === nothing
    odd = FrequencySeries(copy(parent(h)); grid_metadata(h)..., time_samples = 9)
    @test_throws ArgumentError check_same_grid(h, odd)
    @test_throws DimensionMismatch FrequencySeries(copy(parent(h)); grid_metadata(h)..., time_samples = 16)
    @test_throws ArgumentError FrequencySeries(copy(parent(h)); grid_metadata(h)..., time_samples = 1)
    # Arithmetic and full-row slices keep the record; partial bands drop it.
    @test time_samples(h .+ h) == 8
    @test time_samples(h[:, :]) == 8
    @test time_samples(2h) == 8
    @test !(h[2:end, :] isa FrequencySeries)
end

@testset "Series assignment checks the destination grid" begin
    x = TimeSeries(reshape(collect(1.0:16.0), 8, 2); epoch = 100, sampling_interval = 0.25, channels = (:H1, :L1))
    out = copy(x)
    out[:, :] = 2x
    @test parent(out) == 2 .* parent(x)
    shifted = TimeSeries(copy(parent(x)); grid_metadata(x)..., epoch = 101)
    @test_throws ArgumentError (out[:, :] = shifted)
    # Label-indexed assignment compares against the selected channel only.
    single = TimeSeries(ones(8, 1); epoch = 100, sampling_interval = 0.25, channels = (:L1,))
    out[Chan = At([:L1])] = single
    @test all(==(1), parent(out)[:, 2])
    wrong = TimeSeries(ones(8, 1); epoch = 100, sampling_interval = 0.25, channels = (:H1,))
    @test_throws ArgumentError (out[Chan = At([:L1])] = wrong)
    # Validation must not copy the destination block.
    big = TimeSeries(zeros(100_000, 2); epoch = 0, sampling_interval = 1, channels = (:H1, :L1))
    src = copy(big)
    big[:, :] = src
    @test @allocated(big[:, :] = src) < 10_000
end

# A storage type no backend implements: the generic AbstractFFTs fallbacks would recurse or fail obscurely.
struct OpaqueMatrix <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end
Base.size(a::OpaqueMatrix) = size(a.data)
Base.getindex(a::OpaqueMatrix, i::Int, j::Int) = a.data[i, j]

@testset "Missing AbstractFFTs backend" begin
    @test GWSeries.backend_loaded(AbstractFFTs.plan_rfft, Tuple{Matrix{Float64}, Int})
    @test GWSeries.backend_loaded(AbstractFFTs.plan_brfft, Tuple{Matrix{ComplexF64}, Int, Int})
    opaque = OpaqueMatrix(randn(MersenneTwister(2), 8, 1))
    @test !GWSeries.backend_loaded(AbstractFFTs.plan_rfft, Tuple{OpaqueMatrix, Int})
    @test !GWSeries.backend_loaded(AbstractFFTs.plan_brfft, Tuple{OpaqueMatrix, Int, Int})
    @test_throws "no AbstractFFTs backend" GWSeries.require_backend(opaque, false)
    @test_throws "no AbstractFFTs backend" GWSeries.require_backend(opaque, true)
    # rfft reads the wrapped storage directly; the plans and irfft copy into ordinary arrays first.
    x = TimeSeries(opaque; epoch = 0, sampling_interval = 0.25, channels = (:H1,))
    @test_throws "no AbstractFFTs backend" rfft(x)
end
