using Test
import SurrogateWaveforms
import SurrogatePyCallWaveforms
import AbstractWaveforms

@testset "SurrogateWaveforms against PyCall oracle" begin
    params = (
        total_mass = 60.0, luminosity_distance = 500.0, mass_ratio = 2.0,
        spin_1 = (0.2, -0.1, 0.3), spin_2 = (-0.1, 0.2, -0.2),
        inclination = 0.7, phi_ref = 0.3,
        init_orbphase = 0.2, init_quat = (cos(0.1), sin(0.1), 0.0, 0.0),
    )
    for (model, ell_max) in ((:NRSur7dq4, 4), (:NRSur7dq4v2, 5))
        @testset "$model" begin
            surrogate = SurrogateWaveforms.load_surrogate(model)
            native = SurrogateWaveforms.PrecessingGenerator(surrogate; ell_max, f_ref = 25.0)
            oracle = SurrogatePyCallWaveforms.load_generator(SurrogateWaveforms.data_path(model); model, ell_max, f_ref = 25.0)
            labels = Tuple((l, m) for l in 2:ell_max for m in (-l):l)
            for grid in ((;), (; dt = 1 / 2048), (; times = -0.1:(1 / 2048):0.01))
                @testset "$(keys(grid))" begin
                    expected = AbstractWaveforms.waveform_modes_td(oracle, params; grid...)
                    actual = AbstractWaveforms.waveform_modes_td(native, params; grid...)
                    @test actual.t ≈ expected.t
                    @test actual.mode_indices == expected.mode_indices == labels
                    @test size(actual.modes) == size(expected.modes)
                    @test all(isfinite, actual.modes)
                    # Upstream's C++ spline mis-evaluates isolated points in the last
                    # coorbital interval. Use the same exclusion as the fixture tests.
                    tscale = params.total_mass * SurrogateWaveforms.MSUN_IN_SECONDS
                    keep = expected.t .<= surrogate.coorbital.t[end - 1] * tscale
                    @test any(keep)
                    @test maximum(abs, actual.modes[keep, :] - expected.modes[keep, :]) < 1.0e-11 * maximum(abs, expected.modes)
                    isempty(grid) && continue
                    expected_hp = AbstractWaveforms.hphc_td(oracle, params; grid...)
                    actual_hp = AbstractWaveforms.hphc_td(native, params; grid...)
                    @test collect(times(actual_hp)) ≈ collect(times(expected_hp))
                    @test maximum(abs, actual_hp[:, :plus][keep] - expected_hp[:, :plus][keep]) < 1.0e-11 * maximum(abs, expected_hp[:, :plus])
                    @test maximum(abs, actual_hp[:, :cross][keep] - expected_hp[:, :cross][keep]) < 1.0e-11 * maximum(abs, expected_hp[:, :cross])
                end
            end
        end
    end
end
