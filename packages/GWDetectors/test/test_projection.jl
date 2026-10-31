using TestItems

@testsnippet Setup_projection begin
    using StaticArrays

    include("prototypes.jl")
    # A tiny shifted-axis vector verifies that equal lengths do not suffice.
    struct ShiftedVector{T} <: AbstractVector{T}
        data::Vector{T}
    end
    Base.size(v::ShiftedVector) = size(v.data)
    Base.axes(v::ShiftedVector) = (0:(length(v.data) - 1),)
    Base.getindex(v::ShiftedVector, i::Int) = v.data[i + 1]
end

@testitem "projection" setup = [Common, Setup_projection] begin
    d = detectors.H1
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    start_time = sky.geocent_time - 0.125
    for f in (range(0, 100; length = 11), range(0.1; step = 0.3, length = 11), 12.5:23.0:104.5, range(0.0; step = 1, length = 0))
        hp = complex.(1 .+ f / 100, f / 200)
        hc = complex.(f / 300, 2 .- f / 100)
        p = (; plus = hp, cross = hc)
        a = antenna_pattern(sky, d)
        dt = (sky.geocent_time - start_time) + time_delay_from_geocenter(sky, d)
        expected = (a.plus .* hp .+ a.cross .* hc) .* exp.(-2π * im .* f .* dt)
        h = project_fd(p, sky, d; frequencies = f, epoch = start_time)
        @test h ≈ expected atol = 1.0e-13
        @test h !== hp
        @test eltype(h) <: Complex
        @test project_fd((; plus = 2hp, cross = 2hc), sky, d; frequencies = f, epoch = start_time) ≈ 2h
        isempty(f) || first(f) != 0 || @test h[1] ≈ a.plus * hp[1] + a.cross * hc[1]
        ref = sky.geocent_time + 1000
        refsky = merge(sky, (; geocent_time = ref))
        r = antenna_pattern(refsky, d)
        expected_ref = (r.plus .* hp .+ r.cross .* hc) .* exp.(-2π * im .* f .* dt)
        @test project_fd(p, sky, d; frequencies = f, epoch = start_time, reference_time = ref) ≈ expected_ref atol = 1.0e-13
    end
    domain = frequency_prototype(; frequencies = 0:20:20, epoch = start_time, channels = (:H1,))
    @test_throws DimensionMismatch check_data(domain, ones(1))
    @test_throws ArgumentError check_data(domain, ShiftedVector(ones(2)))
    @test_throws ArgumentError check_data(domain, 1)
end
