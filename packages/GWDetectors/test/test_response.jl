using TestItems

@testsnippet Setup_response begin
    using StaticArrays
end

@testitem "response and timing" setup = [Common, Setup_response] begin
    d = Detector((299792458, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1))
    # Source overhead: LAL basis X=-y, Y=-x at psi=0.
    @test antenna_pattern(0, π / 2, 0, 0, d).plus ≈ -1
    @test abs(antenna_pattern(0, π / 2, 0, 0, d).cross) < 1.0e-15
    @test time_delay_from_geocenter(0, 0, 0, d) == -1
    @test time_delay_from_geocenter(π, 0, 0, d) == 1
    swapped = Detector(d.vertex, d.yarm, d.xarm, d.arm_lengths)
    for ra in (0.0, 1.1, 6.2), dec in (-1.4, 0.1, 1.5), psi in (-0.5, 0.0, 2.0)
        f = antenna_pattern(ra, dec, psi, 2.5, d)
        period = antenna_pattern(ra, dec, psi + π, 2.5, d)
        exchange = antenna_pattern(ra, dec, psi, 2.5, swapped)
        @test collect(values(f)) ≈ collect(values(period)) atol = 2.0e-15
        @test collect(values(f)) ≈ -collect(values(exchange)) atol = 2.0e-15
        @test sum(abs2, values(f)) <= 1 + 1.0e-14
        sky = (; ra, dec, geocent_time = 1126259462.25)
        a, b = detectors.H1, detectors.L1
        @test time_delay(sky, a, b) ≈ -time_delay(sky, b, a)
        @test time_delay(sky, a, b) ≈ time_delay_from_geocenter(sky, a) - time_delay_from_geocenter(sky, b) atol = 1.0e-17
        @test abs(time_delay(sky, a, b)) <= light_travel_time(a, b)
        @test light_travel_time(a, b) == light_travel_time(b, a)
        @test light_travel_time(a, a) == 0
    end
    for t in (Float32(1126259462), Float16(100))
        @test greenwich_mean_sidereal_time(t) == greenwich_mean_sidereal_time(Float64(t))
    end
    t = 1126259462.25
    rate = 2π / 86164.0905
    @test greenwich_mean_sidereal_time(t) > 2π
    @test (greenwich_mean_sidereal_time(t + 0.25) - greenwich_mean_sidereal_time(t)) / 0.25 ≈ rate rtol = 1.0e-6
    for leap in GWDetectors.LEAP_SECONDS
        @test greenwich_mean_sidereal_time(leap) == greenwich_mean_sidereal_time(leap + 1)
        @test greenwich_mean_sidereal_time(leap - 1) < greenwich_mean_sidereal_time(leap)
    end
end
