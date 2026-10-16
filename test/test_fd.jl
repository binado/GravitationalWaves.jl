@testset "fd" begin
    LW = LALWaveforms
    f = range(0, 1024; step = 0.25)
    p = (mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3)

    @testset "struct layout" begin
        @test sizeof(LW.LIGOTimeGPS) == 8
        @test sizeof(LW.LALUnit) == 30
        @test fieldoffset(LW.COMPLEX16Sequence, 2) == 8
        @test sizeof(LW.COMPLEX16FrequencySeries) == 128
        offs = ntuple(i -> fieldoffset(LW.COMPLEX16FrequencySeries, i), 6)
        @test offs == (0, 64, 72, 80, 88, 120)
    end

    @testset "shape and masking" begin
        s = WaveformSettings("IMRPhenomD"; minimum_frequency = 30.0, maximum_frequency = 500.0)
        h = hphc_fd(f, p, s)
        @test h.plus isa Vector{ComplexF64} && h.cross isa Vector{ComplexF64}
        @test length(h.plus) == length(h.cross) == length(f)
        out = (f .< 30.0) .| (f .> 500.0)
        @test all(iszero, h.plus[out]) && all(iszero, h.cross[out])
        @test any(!iszero, h.plus[.!out])
    end

    @testset "default maximum frequency is the end of the grid" begin
        h = hphc_fd(f, p, WaveformSettings(:IMRPhenomD))
        h2 = hphc_fd(f, p, WaveformSettings(:IMRPhenomD; maximum_frequency = last(f)))
        @test h.plus == h2.plus && h.cross == h2.cross
    end

    @testset "errors" begin
        @test_throws LALError WaveformSettings("NotAnApproximant")
        s = WaveformSettings("IMRPhenomD")
        @test_throws ArgumentError hphc_fd(f, merge(p, (spin = 1.0,)), s)
        @test_throws ArgumentError hphc_fd(f, (mass_1 = 1.0,), s)
        @test_throws ArgumentError hphc_fd(range(1, 1024; step = 0.25), p, s)
        @test_throws ArgumentError hphc_fd(range(0, 1024; step = -0.25), p, s)
        @test_throws LALError hphc_fd(f, merge(p, (mass_2 = -1.0,)), s)
        # A failure must not poison the next call.
        @test any(!iszero, hphc_fd(f, p, s).plus)
    end

    @testset "show" begin
        @test occursin("IMRPhenomD", sprint(show, WaveformSettings("IMRPhenomD")))
    end
end
