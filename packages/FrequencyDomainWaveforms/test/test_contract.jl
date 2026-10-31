using TestItems

@testsnippet Setup_contract begin
    using AbstractWaveformsTests
end

@testitem "generator contract" setup = [Common, Setup_contract] begin
    params = (
        mass_1 = 36.0, mass_2 = 29.0, luminosity_distance = 500.0, theta_jn = 0.4, phase = 1.3,
    )
    f = range(0, 512; step = 0.25)
    for g in (
            TaylorF2(; minimum_frequency = 30.0, maximum_frequency = 400.0),
            IMRPhenomD(; minimum_frequency = 30.0, maximum_frequency = 400.0),
            IMRPhenomXAS(; minimum_frequency = 30.0, maximum_frequency = 400.0),
        )
        test_generator_contract(g; params, f)
    end
end
