@testset "Domain preparation and noise compatibility" begin
    curve = TabulatedPSD([0, 64], [2, 2])
    for n in (127, 128)
        td = TimeDomain(; epoch = 1126259460.0, nsamples = n, sampling_frequency = 128, channels = (:H1,))
        model = GaussianNoise(curve, td)
        @test check_compatible(td, model) === nothing
        @test check_compatible(TimeDomain(; epoch = 0, nsamples = n, sampling_frequency = 128, channels = (:H1,)), model) === nothing
        fd = FrequencyDomain(; frequencies = frequencies(model), epoch = epoch(td), channels = (:H1,))
        @test check_compatible(fd, model) === nothing
        @test check_data(fd, noise_fd(model)) === nothing
        @test check_data(td, noise_td(model)) === nothing
        @test noise_fd!(Xoshiro(42), zeros(ComplexF64, nsamples(fd)), model) == noise_fd(Xoshiro(42), model)
        bad = TimeDomain(; epoch = 0, nsamples = n, sampling_frequency = 64, channels = (:H1,))
        @test_throws ArgumentError check_compatible(bad, model)
        labels = FrequencyDomain(; frequencies = frequencies(model), epoch = 0, channels = (:L1,))
        @test_throws ArgumentError check_compatible(labels, model)
    end
    @test_throws ArgumentError GaussianNoise(curve, TimeDomain(; epoch = 0, nsamples = 128, sampling_frequency = 128, channels = (:H1, :L1)))
    @test_throws ArgumentError GaussianNoise(curve, TimeDomain(; epoch = 0, nsamples = 128, sampling_frequency = big"1e400", channels = (:H1,)))
end
