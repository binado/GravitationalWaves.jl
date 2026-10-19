@testset "Domain projection preconditions" begin
    f = 0:0.25:8
    pol = (; plus = complex.(sin.(f)), cross = complex.(cos.(f)))
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    fd = FrequencyDomain(; frequencies = f, epoch = 1126259460.0, channels = (:H1,))
    @test check_data(fd, pol.plus) === nothing
    h = project_fd(fd, pol, sky, detectors.H1)
    @test h isa Vector{ComplexF64}
    renamed = FrequencyDomain(; frequencies = f, epoch = epoch(fd), channels = (:L1,))
    @test project_fd(renamed, pol, sky, detectors.H1) == h
    @test_throws MethodError project_fd(collect(f), pol, sky, detectors.H1; start_time = epoch(fd))
end
