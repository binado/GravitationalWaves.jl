using FFTW
using LinearAlgebra

@testset "Physical FFTs and reusable plans" begin
    for T in (Float32, Float64), n in (2, 6, 16), nc in (1, 2), viewed in (false, true)
        a = reshape(T.(sin.(1:(n * nc))), n, nc)
        nc == 1 && (a = vec(a))
        viewed && (a = view(a, ntuple(_ -> Colon(), ndims(a))...))
        saved = copy(a)
        names = nc == 1 ? (:H1,) : (:H1, :L1)
        td = TimeDomain(; epoch = 1126259460.125, nsamples = n, sampling_frequency = 32, channels = names)
        fd = frequency_domain(td)
        check_data(td, a)
        h = to_frequencydomain(td, a)
        @test a == saved
        @test eltype(h) == Complex{T}
        @test h ≈ rfft(a, 1) ./ 32
        @test check_same_grid(time_domain(fd), td) === nothing
        h_saved = copy(h)
        x = to_timedomain(fd, h)
        @test h == h_saved
        @test eltype(x) == T
        @test x ≈ a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        forward = plan_to_frequencydomain(td, a)
        inverse = plan_to_timedomain(fd, h)
        @test check_same_grid(output_domain(forward), fd) === nothing
        @test check_same_grid(output_domain(inverse), td) === nothing
        @test transform(forward, a) ≈ h
        @test transform(inverse, h) ≈ a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        desth, destx = similar(h), similar(x)
        @test transform!(desth, forward, a) === desth
        @test desth ≈ h
        @test transform!(destx, inverse, h) === destx
        @test destx ≈ a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        @test a == saved && h == h_saved
        @test transform!(desth, forward, 2 .* a) ≈ 2h
        @test transform!(destx, inverse, 2 .* h) ≈ 2a rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        shifted = TimeDomain(; epoch = 0, nsamples = n, sampling_frequency = 32, channels = names)
        @test to_frequencydomain(shifted, a) == h
    end
    td = TimeDomain(; epoch = 0, nsamples = 16, sampling_frequency = 16, channels = (:H1,))
    fd = frequency_domain(td)
    a = collect(1.0:16.0)
    h = to_frequencydomain(td, a)
    power = abs2.(h)
    @test sum(abs2, a) * sample_spacing(td) ≈ frequency_spacing(fd) * (power[1] + power[end] + 2sum(power[2:(end - 1)]))
    low, vals = crop(fd, h; start = 1, stop = 9)
    @test_throws ArgumentError time_domain(low)
    restored, vals = pad(low, vals; before = 1)
    @test nsamples(time_domain(restored)) == 16
    high, vals = crop(fd, h; start = 0, stop = 5)
    @test sampling_frequency(time_domain(high)) == 8
    @test duration(time_domain(high)) == 1
    restored, vals = pad(high, vals; after = 4)
    @test sampling_frequency(time_domain(restored)) == 16
    dirty = copy(h);dirty[1] += 7im;dirty[end] -= 9im
    @test to_timedomain(fd, dirty) == to_timedomain(fd, h)
    inverse = plan_to_timedomain(fd, dirty)
    @test transform(inverse, dirty) == to_timedomain(fd, h)
    @test dirty[1] == h[1] + 7im && dirty[end] == h[end] - 9im
    realfd = FrequencyDomain(; frequencies = 0:1:4, epoch = 0, channels = (:H1,))
    @test transform(plan_to_timedomain(realfd, ones(Float32, 5)), ones(Float32, 5)) ≈ Float32[8, 0, 0, 0, 0, 0, 0, 0]
    precise = TimeDomain(; epoch = big"1126259460.123456789", nsamples = 4, sampling_frequency = big"32", channels = (:H1,))
    @test epoch(time_domain(frequency_domain(precise))) == epoch(precise)
    @test_throws DimensionMismatch plan_to_frequencydomain(td, ones(8))
    odd = TimeDomain(; epoch = 0, nsamples = 3, sampling_frequency = 4, channels = (:H1,))
    @test_throws ArgumentError plan_to_frequencydomain(odd, ones(3))
    @test_throws ArgumentError plan_to_timedomain(low, h[2:end])
    indexed = view(a, [1, 3, 4, 8])
    @test_throws MethodError to_frequencydomain(precise, indexed)
end

@testset "Nonbinary sampling rate and strain-scale samples" begin
    for T in (Float32, Float64)
        td = TimeDomain(; epoch = 1126259460.125, nsamples = 14, sampling_frequency = 100.1, channels = (:H1, :L1, :V1))
        fd = frequency_domain(td)
        x = reshape(T.(sin.(1:42) .* 1.0e-22), 14, 3)
        h = to_frequencydomain(td, x)
        @test h ≈ rfft(x, 1) .* sample_spacing(td) rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        forward = plan_to_frequencydomain(td, x)
        inverse = plan_to_timedomain(fd, h)
        @test transform(forward, x) ≈ h
        @test transform(inverse, h) ≈ to_timedomain(fd, h)
        @test transform(inverse, h) ≈ x rtol = (T == Float32 ? 1.0e-6 : 1.0e-14)
        @test check_same_grid(output_domain(inverse), td) === nothing
    end
end

@testset "Strided FFT prototypes" begin
    td = TimeDomain(; epoch = 0, nsamples = 14, sampling_frequency = 100.1, channels = (:H1, :L1, :V1))
    fd = frequency_domain(td)
    backing = reshape(sin.(1:84), 28, 3)
    samples = view(backing, 1:2:28, :)
    saved = copy(samples)
    forward = plan_to_frequencydomain(td, samples)
    h = transform(forward, samples)
    @test h ≈ to_frequencydomain(td, samples)
    bins = zeros(ComplexF64, 2nsamples(fd), 3)
    strided = view(bins, 1:2:size(bins, 1), :)
    strided .= h
    inverse = plan_to_timedomain(fd, strided)
    @test transform(inverse, strided) ≈ samples
    @test to_timedomain(fd, strided) ≈ samples
    @test samples == saved && strided == h
end
