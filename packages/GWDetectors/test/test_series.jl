using TestItems

@testsnippet Setup_series begin
    using StaticArrays

    include("prototypes.jl")
end

@testitem "Domain projection preconditions" setup = [Common, Setup_series] begin
    f = 0:0.25:8
    pol = (; plus = complex.(sin.(f)), cross = complex.(cos.(f)))
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    fd = frequency_prototype(; frequencies = f, epoch = 1126259460.0, channels = (:H1,))
    @test check_data(fd, pol.plus) === nothing
    h = project_fd(pol, sky, detectors.H1; frequencies = frequencies(fd), epoch = epoch(fd))
    @test h isa Vector{ComplexF64}
    renamed = frequency_prototype(; frequencies = f, epoch = epoch(fd), channels = (:L1,))
    @test project_fd(pol, sky, detectors.H1; frequencies = frequencies(renamed), epoch = epoch(renamed)) == h
    series = FrequencySeries(hcat(pol.plus, pol.cross); frequency_spacing = step(f), epoch = 0, channels = (:plus, :cross))
    actual = project_fd(series, sky, detectors.H1; epoch = epoch(fd), channel = :H1)
    @test actual isa FrequencySeries
    @test channels(actual) == (:H1,)
    @test parent(actual)[:, 1] ≈ h
    # A source transform with a negative start must retain its phase reference.
    shifted = FrequencySeries(parent(series) .* exp.(2π * im .* f .* -0.125); grid_metadata(series)..., epoch = -0.125)
    @test parent(project_fd(shifted, sky, detectors.H1; epoch = epoch(fd))) ≈ parent(actual)
end
