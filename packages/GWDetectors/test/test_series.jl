using GWSeries

@testset "Series projection delegation" begin
    f = 0:0.25:8
    pol = (; plus = complex.(sin.(f)), cross = complex.(cos.(f)))
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    h = project_fd(FrequencySeries, f, pol, sky, detectors.H1; start_time = 1126259460.0, channels = (:H1,))
    @test parent(h) == project_fd(f, pol, sky, detectors.H1; start_time = epoch(h))
    @test channels(h) == (:H1,)
    @test frequencies(h) == f
    @test epoch(h) == 1126259460.0
    @test_throws DimensionMismatch project_fd(FrequencySeries, f, pol, sky, detectors.H1; start_time = 0, channels = (:H1, :L1))
    @test_throws ArgumentError project_fd(FrequencySeries, -8:0.25:0, pol, sky, detectors.H1; start_time = 0, channels = (:H1,))
    @test_throws ArgumentError project_fd(FrequencySeries, f, pol, sky, detectors.H1; start_time = Inf, channels = (:H1,))
    @test_throws DimensionMismatch project_fd(FrequencySeries, f, (; plus = pol.plus, cross = [0im]), sky, detectors.H1; start_time = 0, channels = (:H1,))
end
