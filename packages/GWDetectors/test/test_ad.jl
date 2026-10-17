include("ad_helpers.jl")

@testset "ForwardDiff" begin
    for loss in (geometry_loss, antenna_loss, timing_loss, projection_loss)
        actual = ForwardDiff.gradient(loss, AD_INPUT)
        expected = map(eachindex(AD_INPUT)) do i
            # Binary steps preserve the epoch difference exactly in Float64.
            step = i in (4, 5) ? 1 / 8192 : 1.0e-5
            values = map((-2, -1, 1, 2)) do k
                shifted = copy(AD_INPUT)
                shifted[i] += k * step
                loss(shifted)
            end
            (values[1] - 8values[2] + 8values[3] - values[4]) / (12step)
        end
        @test all(isfinite, actual)
        # The five-point stencil resolves the high-frequency timing phase.
        @test actual ≈ expected rtol = 2.0e-6 atol = 5.0e-7
        @test actual[6:11] ≈ expected[6:11] rtol = 1.0e-6 atol = 5.0e-7
    end
    @test ForwardDiff.derivative(greenwich_mean_sidereal_time, 1126259462.25) ≈ 2π / 86164.0905 rtol = 1.0e-8
end
