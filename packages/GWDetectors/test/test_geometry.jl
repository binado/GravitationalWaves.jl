@testset "geometry" begin
    @test keys(detectors) == (:H1, :L1, :V1, :K1, :E1, :E2, :E3, :S1, :S2, :R1, :R2, :C1, :C2)
    for (key, d) in pairs(detectors)
        @test d isa Detector
        @test detectors[key] === d
        @test all(>(0), d.arm_lengths)
        D = detector_tensor(d)
        @test D isa SMatrix{3, 3}
        @test all(D[i, j] == D[j, i] for i in 1:3, j in 1:3)
        @test abs(sum(D[i, i] for i in 1:3)) < 1.0e-6
    end
    @test detectors.H1.vertex == SVector(-2161414.92636, -3834695.17889, 4600350.22664)
    @test detectors.H1.arm_lengths == SVector(3995.084, 3995.044)
    @test_throws UndefKeywordError Detector()
    args = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1000, 1000))
    d = Detector(args...)
    @test detector_tensor(d) == @SMatrix [0.5 0 0; 0 -0.5 0; 0 0 0]
    @test Detector(map(x -> Float32.(x), args)...).vertex isa SVector{3, Float32}
    @test Detector(map(x -> BigFloat.(x), args)...).vertex isa SVector{3, BigFloat}
    @test_throws DimensionMismatch Detector((0, 0), args[2:end]...; check = false)
    @test_throws DimensionMismatch Detector(args[1:3]..., (1,); check = false)
    for (vertex, u, v, lengths) in (
            ((NaN, 0, 0), args[2:end]...), (args[1], (2, 0, 0), args[3:4]...),
            (args[1:2]..., args[2], args[4]), (args[1:3]..., (0, 1)),
            (args[1], (1 + im, 0, 0), args[3:4]...),
        )
        @test_throws ArgumentError Detector(vertex, u, v, lengths)
    end
    @test Detector(args[1], (2, 0, 0), args[3:4]...; check = false).xarm == SVector(2, 0, 0)
    equator = Detector(;
        latitude = 0.0, longitude = 0.0,
        xarm_azimuth = 0.0, yarm_azimuth = π / 3, arm_lengths = (1000, 1000)
    )
    @test equator.vertex == SVector(6378137.0, 0, 0)
    @test equator.xarm == SVector(0, 0, 1)
    @test equator.yarm[2] ≈ sin(π / 3)
    @test dot(equator.xarm, equator.yarm) ≈ 0.5
    tilted = Detector(;
        latitude = 0.0, longitude = 0.0,
        xarm_azimuth = 0.0, yarm_azimuth = π / 2, xarm_tilt = 0.1,
        yarm_tilt = -0.2, elevation = 20, arm_lengths = (3000, 4000)
    )
    @test tilted.vertex == SVector(6378157, 0, 0)
    @test tilted.xarm[1] ≈ sin(0.1)
    @test tilted.yarm[1] ≈ sin(-0.2)
    @test abs(sum(detector_tensor(tilted)[i, i] for i in 1:3)) < 1.0e-15
    for kwargs in ((; latitude = 2.0), (; xarm_tilt = 2.0), (; longitude = Inf), (; elevation = Inf))
        base = (; latitude = 0.0, longitude = 0.0, xarm_azimuth = 0.0, yarm_azimuth = 1.0, arm_lengths = (1, 1))
        @test_throws ArgumentError Detector(; merge(base, kwargs)...)
    end
    @test_throws UndefKeywordError Detector(; latitude = 0.0)
    @test_throws UndefKeywordError Detector(; elevation = Inf)
end
