using TestItems

@testsnippet Setup_quality begin
    using Aqua
    using JET
end

@testitem "Aqua" setup = [Common, Setup_quality] begin
    Aqua.test_all(SurrogateWaveforms)
end

@testitem "JET" setup = [Common, Setup_quality] begin
    # JET's analysis is tied to the Julia version it was built for.
    if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
        @testset "JET" begin
            JET.test_package(SurrogateWaveforms; target_modules = (SurrogateWaveforms,))
        end

        # The kernels must be free of runtime dispatch for the AD backends.
        @testset "JET optimisation" begin
            s = SurrogateWaveforms.NaturalCubicSpline([0.0, 1.0, 2.5, 4.0], [0.0, 1.0, 0.5, 2.0])
            JET.@test_opt target_modules = (SurrogateWaveforms,) SurrogateWaveforms.spline_eval(s, 1.3)
            JET.@test_opt target_modules = (SurrogateWaveforms,) SurrogateWaveforms.wigner_D_element(3, 1, -2, 0.6 + 0.1im, 0.3 - 0.2im)
            JET.@test_opt target_modules = (SurrogateWaveforms,) SurrogateWaveforms.rescale_to_norm((0.1, 0.2, 0.3), 0.5)
        end
    end
end

@testitem "type stability" setup = [Common, Setup_quality] begin
    s = SurrogateWaveforms.NaturalCubicSpline([0.0, 1.0, 2.5, 4.0], [0.0, 1.0, 0.5, 2.0])
    @test (@inferred SurrogateWaveforms.spline_eval(s, 1.3)) isa Float64
    @test (@inferred SurrogateWaveforms.wigner_D_element(3, 1, -2, 0.6 + 0.1im, 0.3 - 0.2im)) isa ComplexF64
    @test (@inferred SurrogateWaveforms.rescale_to_norm((0.1, 0.2, 0.3), 0.5)) isa NTuple{3, Float64}
end
