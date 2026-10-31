using TestItems

@testsnippet Setup_sensitivity begin
    using StaticArrays
end

@testitem "detector sensitivities" setup = [Common, Setup_sensitivity] begin
    @test keys(sensitivities) == keys(detectors)
    @test sensitivities.H1 == sensitivities.L1 == :AplusDesign
    @test sensitivities.V1 == :avirgo_O5low_NEW
    @test sensitivities.K1 == :kagra_128Mpc
    for name in (:E1, :E2, :E3)
        @test sensitivities[name] == :ET_COBA_10km
        @test all(==(10000), detectors[name].arm_lengths)
    end
    for name in (:S1, :S2, :R1, :R2)
        @test sensitivities[name] == :ET_COBA_15km
        @test all(==(15000), detectors[name].arm_lengths)
    end
    for name in (:C1, :C2)
        @test sensitivities[name] == :cosmic_explorer_strain
        @test all(==(40000), detectors[name].arm_lengths)
    end
end
