using DataInterpolations
using FiniteDifferences
using ForwardDiff

@testset "NaturalCubicSpline" begin
    t = [0.0, 0.3, 0.9, 1.4, 2.2, 3.0, 3.1, 4.5]
    xs = range(0.0, 4.5; length = 53)

    @testset "matches DataInterpolations.CubicSpline" begin
        u = @. sin(2t) + 0.1t^2
        ref = DataInterpolations.CubicSpline(u, t)
        s = GWSurrogates.NaturalCubicSpline(t, u)
        @test maximum(abs(GWSurrogates.spline_eval(s, x) - ref(x)) for x in xs) < 1.0e-13

        uc = complex.(u, cos.(3t))
        refc = DataInterpolations.CubicSpline(uc, t)
        sc = GWSurrogates.NaturalCubicSpline(t, uc)
        @test maximum(abs(GWSurrogates.spline_eval(sc, x) - refc(x)) for x in xs) < 1.0e-13
    end

    @testset "interpolates the knots and is natural" begin
        u = @. exp(-t)
        s = GWSurrogates.NaturalCubicSpline(t, u)
        @test [GWSurrogates.spline_eval(s, ti) for ti in t] ≈ u atol = 1.0e-14
        @test s.m[1] == 0 && s.m[end] == 0
    end

    @testset "reproduces a line and rejects bad input" begin
        s = GWSurrogates.NaturalCubicSpline(t, 2 .* t .+ 1)
        @test all(x -> GWSurrogates.spline_eval(s, x) ≈ 2x + 1, xs)
        @test_throws DimensionMismatch GWSurrogates.NaturalCubicSpline(t, zeros(3))
        @test_throws ArgumentError GWSurrogates.NaturalCubicSpline([0.0, 1.0], [0.0, 1.0])
        @test_throws ArgumentError GWSurrogates.NaturalCubicSpline([0.0, 1.0, 1.0], zeros(3))
    end

    @testset "ForwardDiff agrees with finite differences" begin
        u0 = @. sin(2t) + 0.1t^2
        fdm = central_fdm(5, 1)
        x0 = 1.234
        # Gradient with respect to the data values.
        f_u(u) = GWSurrogates.spline_eval(GWSurrogates.NaturalCubicSpline(t, u), x0)
        @test ForwardDiff.gradient(f_u, u0) ≈ grad(fdm, f_u, u0)[1] atol = 1.0e-8
        # Derivative with respect to the evaluation point, away from knots.
        s = GWSurrogates.NaturalCubicSpline(t, u0)
        f_x(x) = GWSurrogates.spline_eval(s, x)
        for x in (0.1, 1.234, 2.7, 4.0)
            @test ForwardDiff.derivative(f_x, x) ≈ fdm(f_x, x) atol = 1.0e-8
        end
        # Gradient with respect to the knot positions.
        f_t(t) = GWSurrogates.spline_eval(GWSurrogates.NaturalCubicSpline(t, u0), x0)
        @test ForwardDiff.gradient(f_t, t) ≈ grad(fdm, f_t, t)[1] atol = 1.0e-7
    end
end
