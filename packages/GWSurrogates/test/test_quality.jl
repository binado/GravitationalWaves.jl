using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(GWSurrogates)
end

# JET's analysis is tied to the Julia version it was built for.
if VERSION.major == 1 && VERSION.minor == 12 && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(GWSurrogates; target_modules = (GWSurrogates,))
    end

    # The kernels must be free of runtime dispatch for the AD backends.
    @testset "JET optimisation" begin
        s = GWSurrogates.NaturalCubicSpline([0.0, 1.0, 2.5, 4.0], [0.0, 1.0, 0.5, 2.0])
        JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.spline_eval(s, 1.3)
        JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.wigner_D_element(3, 1, -2, 0.6 + 0.1im, 0.3 - 0.2im)
        JET.@test_opt target_modules = (GWSurrogates,) GWSurrogates.rescale_to_norm((0.1, 0.2, 0.3), 0.5)
    end
end

@testset "type stability" begin
    s = GWSurrogates.NaturalCubicSpline([0.0, 1.0, 2.5, 4.0], [0.0, 1.0, 0.5, 2.0])
    @test (@inferred GWSurrogates.spline_eval(s, 1.3)) isa Float64
    @test (@inferred GWSurrogates.wigner_D_element(3, 1, -2, 0.6 + 0.1im, 0.3 - 0.2im)) isa ComplexF64
    @test (@inferred GWSurrogates.rescale_to_norm((0.1, 0.2, 0.3), 0.5)) isa NTuple{3, Float64}
end
