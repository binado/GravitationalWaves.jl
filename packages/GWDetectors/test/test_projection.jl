# A tiny shifted-axis vector verifies that equal lengths do not suffice.
struct ShiftedVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(v::ShiftedVector) = size(v.data)
Base.axes(v::ShiftedVector) = (0:(length(v.data) - 1),)
Base.getindex(v::ShiftedVector, i::Int) = v.data[i + 1]

@testset "projection" begin
    d = Detector("H1")
    sky = (; ra = 1.2, dec = -0.4, psi = 0.7, geocent_time = 1126259462.25)
    start_time = sky.geocent_time - 0.125
    for f in (range(0, 100; length = 11), [0.0, 12.3, 1.2, 90.0], Float64[])
        hp = complex.(1 .+ f / 100, f / 200)
        hc = complex.(f / 300, 2 .- f / 100)
        p = (; plus = hp, cross = hc)
        a = antenna_pattern(sky, d)
        dt = (sky.geocent_time - start_time) + time_delay_from_geocenter(sky, d)
        expected = (a.plus .* hp .+ a.cross .* hc) .* exp.(-2π * im .* f .* dt)
        h = project_fd(f, p, sky, d; start_time)
        @test h ≈ expected atol = 1.0e-13
        @test h !== hp
        @test eltype(h) <: Complex
        @test project_fd(f, (; plus = 2hp, cross = 2hc), sky, d; start_time) ≈ 2h
        @test h == project_fd(f, p, sky, d; start_time, check = false)
        isempty(f) || @test h[1] ≈ a.plus * hp[1] + a.cross * hc[1]
        ref = sky.geocent_time + 1000
        refsky = merge(sky, (; geocent_time = ref))
        r = antenna_pattern(refsky, d)
        expected_ref = (r.plus .* hp .+ r.cross .* hc) .* exp.(-2π * im .* f .* dt)
        @test project_fd(f, p, sky, d; start_time, reference_time = ref) ≈ expected_ref atol = 1.0e-13
    end
    f = [0.0, 20.0]
    p = (; plus = ones(ComplexF64, 2), cross = zeros(ComplexF64, 2))
    @test_throws DimensionMismatch project_fd(f, (; plus = ones(1), cross = ones(2)), sky, d; start_time, check = false)
    @test_throws DimensionMismatch project_fd(f, (; plus = ShiftedVector(ones(2)), cross = ones(2)), sky, d; start_time, check = false)
    @test_throws ArgumentError project_fd(f, (; plus = 1, cross = 1), sky, d; start_time, check = false)
    for bad in ([-1.0, 20], [NaN, 20], [Inf, 20])
        @test_throws ArgumentError project_fd(bad, p, sky, d; start_time)
    end
    @test_throws ArgumentError project_fd(f, p, sky, d; start_time = NaN)
    @test_throws ArgumentError project_fd(f, p, sky, d; start_time, reference_time = Inf)
    @test_throws ArgumentError project_fd(f, (; plus = [NaN, 1], cross = ones(2)), sky, d; start_time)
end
