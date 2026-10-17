@testset "generator contract" begin
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
    )
    f = range(0, 512; step = 0.25)
    # IMRPhenomD is native FD, IMRPhenomXPHM is precessing with a mode array, TaylorT4 goes
    # through the time-domain path.
    for name in ("IMRPhenomD", "IMRPhenomXPHM", "TaylorT4")
        test_generator_contract(
            LALWaveformGenerator(name; minimum_frequency = 30.0, maximum_frequency = 400.0);
            params, f,
        )
    end
end
