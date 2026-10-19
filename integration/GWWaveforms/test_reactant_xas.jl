using Reactant, Enzyme, ForwardDiff
Reactant.set_default_backend("cpu")

@testset "IMRPhenomXAS Reactant CPU" begin
    grid = range(0, 512; step = 0.5)
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 450.0, theta_jn = 0.7, phase = 1.1,
        a_1 = 0.3, a_2 = 0.2, tilt_1 = 0.0, tilt_2 = Float64(π),
    )
    draws = (
        params, merge(params, (mass_1 = 30.0, mass_2 = 30.0)),
        merge(params, (mass_1 = 12.0, mass_2 = 31.0, a_1 = 0.8, a_2 = 0.6)),
        merge(params, (mass_1 = 180.0, mass_2 = 120.0)),
    )
    @testset "one executable, moving regions and cutoff" begin
        for options in ((;), (; minimum_frequency = 20.1, maximum_frequency = 200.3, reference_frequency = 0.0))
            g = IMRPhenomXAS(; options...)
            tracked = Reactant.to_rarray(params; track_numbers = true)
            executable = Reactant.@compile hphc_fd(grid, tracked, g)
            for p in draws
                validate_parameters(grid, p, g)
                actual = executable(grid, Reactant.to_rarray(p; track_numbers = true), g)
                expected = hphc_fd(grid, p, g)
                @test maximum(abs, Array(actual.plus) - expected.plus) <= 1.0e-10 * maximum(abs, expected.plus)
                @test maximum(abs, Array(actual.cross) - expected.cross) <= 1.0e-10 * maximum(abs, expected.cross)
                @test iszero.(Array(actual.plus)) == iszero.(expected.plus)
            end
        end
    end
    @testset "compiled waveform gradients" begin
        g = IMRPhenomXAS(minimum_frequency = 25.0)
        template = 1.0e22 .* hphc_fd(AD_GRID, ad_params(AD_X0 .* [1.01, 1, 1, 1, 1, 1, 1]), g).cross
        loss = ad_loss(g, template)
        objective(x) = Reactant.@allowscalar loss(x)
        executable = Reactant.@compile Enzyme.gradient(Enzyme.Reverse, objective, Reactant.to_rarray(AD_X0))
        for x in (AD_X0, [130.0, 0.8, 0.3, 0.2, 450.0, 0.7, 1.1], AD_X0 .* [1.05, 0.9, 1, 1, 1.1, 0.8, 0.5])
            expected = ForwardDiff.gradient(loss, x)
            actual = Array(only(executable(Enzyme.Reverse, objective, Reactant.to_rarray(x))))
            @test all(isfinite, actual)
            @test actual ≈ expected rtol = 1.0e-8
        end
    end
end
