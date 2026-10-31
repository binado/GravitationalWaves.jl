using TestItems

@testitem "bundled PSDs" setup = [Common] begin
    @test length(available_psds()) == 12
    @test length(unique(available_psds())) == 12
    for name in available_psds()
        @test psd_path(String(name)) == psd_path(name)
        curve = TabulatedPSD(psd_path(name), ' ')
        @test all(>(0), curve.values)
        i = length(curve.frequencies) ÷ 2
        @test psd(curve.frequencies[i], curve) ≈ curve.values[i]
        @test_throws DomainError psd(last(curve.frequencies) + 1, curve)
    end
    @test_throws ArgumentError psd_path(:unknown)
    @test_throws ArgumentError psd_path("../geometry.toml")
end
