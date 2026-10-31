using TestItems

@testsnippet Setup_quality begin
    using StaticArrays

    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(GWDetectors)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(GWDetectors; target_modules = (GWDetectors,))
            @test (@inferred antenna_pattern(1.2, -0.4, 0.7, 2.0, detectors.H1)).plus isa Float64
            @test (@inferred greenwich_mean_sidereal_time(1126259462.25)) isa Float64
        end
    end
end
