using Test
using GWTimeFrequencyDomains
using Random
using FFTW: rfft, irfft
using GWNoise

@testset "GWNoise" begin
    @testset "Spectral density tables" begin
        curve = read_psd(IOBuffer("# frequency PSD\n\n0 2\n2 6 # comment\n4 10\n"))
        @test psd([0, 1, 4], curve) == [2, 4, 10]
        @test asd(1, curve) == 2
        @test psd([-1, 5], curve; outside = :zero) == [0, 0]
        @test_throws DomainError psd(5, curve)
        @test_throws ArgumentError TabulatedPSD([1, 1], [2, 3])
        @test_throws ArgumentError TabulatedPSD([1, 2], [2, -3])
        mktemp() do path, io
            write(io, "# ASD\n1,2\n3,4\n")
            flush(io)
            curve = read_asd(path)
            @test psd([1, 2, 3], curve) == [4, 10, 16]
            @test asd([1, 3], curve) == [2, 4]
        end
        powerlaw = TabulatedPSD([1, 10, 100], [1, 100, 10000]; interpolation = :loglog)
        @test psd([2, 20, 50], powerlaw) ≈ [4, 400, 2500]
        @test_throws ArgumentError TabulatedPSD([0, 1], [1, 2]; interpolation = :loglog)
    end

    @testset "Reading files" begin
        loglog = read_psd(IOBuffer("1 1\n10 100\n100 10000\n"); interpolation = :loglog)
        @test psd(2, loglog) ≈ 4
        csv = read_psd(IOBuffer("# frequency, PSD\n1, 2\n3 ,4\n"))
        @test psd([1, 3], csv) == [2, 4]
        tabs = read_psd(IOBuffer("1\t2\n3\t4\n"))
        @test psd([1, 3], tabs) == [2, 4]
        @test_throws ArgumentError read_psd(IOBuffer("1 2 3\n4 5 6\n"))
        @test_throws ArgumentError read_psd(IOBuffer(""))
        @test_throws ArgumentError read_psd(IOBuffer("# only a comment\n\n"))
        @test_throws ErrorException read_psd(IOBuffer("1 2\n3 abc\n"))
        @test_throws ArgumentError read_psd(IOBuffer("1 2\n3 NaN\n"))
        @test_throws ArgumentError read_psd(IOBuffer("1 2\n3 -4\n"))
        @test_throws ArgumentError read_asd(IOBuffer("1 Inf\n3 4\n"))
    end

    @testset "Sampling and Fourier convention" begin
        curve = TabulatedPSD([0, 64], [2, 2])
        for n in (127, 128)
            model = GaussianNoise(curve, TimeDomain(; epoch = 0.0, nsamples = n, sampling_frequency = 128, channels = (:H1,)))
            @test length(frequencies(model)) == n ÷ 2 + 1
            @test step(frequencies(model)) == 128 / n
            fd = noise_fd(Xoshiro(42), model)
            td = noise_td(Xoshiro(42), model)
            @test fd[1] == 0
            @test iseven(n) ? fd[end] == 0 : fd[end] != 0
            @test length(td) == n
            @test rfft(td) / 128 ≈ fd
            @test noise_fd!(Xoshiro(42), similar(fd), model) == fd
            @test length(noise_fd(model)) == length(fd)
            @test length(noise_td(model)) == n
            @test noise_fd!(similar(fd), model) isa Vector{ComplexF64}
        end
        model = GaussianNoise(curve, TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)); minimum_frequency = 10.5, maximum_frequency = 30)
        expected = [10.5 <= x <= 30 ? 2.0 : 0.0 for x in frequencies(model)]
        @test psd(model) == expected
        @test all(iszero, noise_fd(model)[expected .== 0])
        @test_throws ArgumentError GaussianNoise(curve, TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)); maximum_frequency = 65)
        @test_throws DimensionMismatch check_data(FrequencyDomain(; frequencies = frequencies(model), epoch = 0, channels = (:H1,)), zeros(ComplexF64, 3))
        @test_throws DimensionMismatch check_data(TimeDomain(; epoch = 0, nsamples = nsamples(model.time_domain), sampling_frequency = sampling_frequency(model.time_domain), channels = (:H1,)), zeros(3))
        @test_throws MethodError noise_fd!(zeros(length(frequencies(model))), model)
        for kwargs in ((sampling_frequency = 0, nsamples = 128), (sampling_frequency = 128, nsamples = 1))
            @test_throws ArgumentError GaussianNoise(curve, TimeDomain(; epoch = 0, channels = (:H1,), kwargs...))
        end
        @test_throws ArgumentError GaussianNoise(TabulatedPSD([20, 100], [2, 4]), TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)); minimum_frequency = 10, maximum_frequency = 50)
        limited = GaussianNoise(TabulatedPSD([20, 100], [2, 4]), TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)))
        @test all(iszero, psd(limited)[frequencies(limited) .< 20])
        @test psd(limited)[21] == 2
        rng = Xoshiro(42)
        @test noise_fd(rng, model) != noise_fd(rng, model)
    end

    @testset "Gaussian normalization" begin
        model = GaussianNoise(TabulatedPSD([0, 64], [2, 2]), TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)))
        rng = Xoshiro(123)
        count = 8000
        draws = reduce(hcat, (noise_fd(rng, model) for _ in 1:count))
        re, im = real.(draws[11, :]), imag.(draws[11, :])
        # S/(4Δf) = 1/2 per quadrature; independent bins and quadratures.
        @test abs(sum(re) / count) < 0.04
        @test sum(abs2, re) / count ≈ 0.5 rtol = 0.06
        @test sum(abs2, im) / count ≈ 0.5 rtol = 0.06
        @test abs(sum(re .* im) / count) < 0.04
        @test abs(sum(conj.(draws[11, :]) .* draws[12, :]) / count) < 0.05
        colored = GaussianNoise(TabulatedPSD([0, 64], [1, 129]), TimeDomain(; epoch = 0.0, nsamples = 128, sampling_frequency = 128, channels = (:H1,)))
        colored_draws = reduce(hcat, (noise_fd(rng, colored) for _ in 1:count))
        for bin in (11, 41)
            @test sum(abs2, colored_draws[bin, :]) / count ≈ psd(colored)[bin] / 2 rtol = 0.06
        end
        td = irfft(draws, nsamples(model.time_domain), 1) .* sampling_frequency(model.time_domain)
        # Each active interior bin contributes S Δf to the TD variance.
        @test sum(abs2, td) / length(td) ≈ sum(psd(model)) * step(frequencies(model)) rtol = 0.02
    end

    include("test_domains.jl")
    include("test_quality.jl")
end
