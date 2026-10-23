using GWSeries
include("../cases.jl")
include("../oracle_checks.jl")
using Test
import SurrogateWaveforms
import SurrogatePyCallWaveforms
import AbstractWaveforms

@testset "SurrogateWaveforms against PyCall oracle" begin
    for (model, ell_max) in ((:NRSur7dq4, 4), (:NRSur7dq4v2, 5)), id in WaveformCases.SURROGATE_IDS
        params = WaveformCases.surrogate_parameters(WaveformCases.CASES[id])
        @testset "$model / $id" begin
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
                    @test all(isfinite, expected.modes)
                    # Upstream's C++ spline mis-evaluates isolated points in the last
                    # coorbital interval. Use the same exclusion as the fixture tests.
                    tscale = params.total_mass * SurrogateWaveforms.MSUN_IN_SECONDS
                    keep = expected.t .<= surrogate.coorbital.t[end - 1] * tscale
                    @test any(keep)
                    check_oracle(actual.modes[keep, :], expected.modes[keep, :]; rtol = 1.0e-11)
                    isempty(grid) && continue
                    expected_hp = AbstractWaveforms.hphc_td(oracle, params; grid...)
                    actual_hp = AbstractWaveforms.hphc_td(native, params; grid...)
                    @test collect(times(actual_hp)) ≈ collect(times(expected_hp))
                    @test all(isfinite, actual_hp) && all(isfinite, expected_hp)
                    @test channels(actual_hp) == channels(expected_hp) == (:plus, :cross)
                    for channel in (:plus, :cross)
                        check_oracle(actual_hp[:, channel][keep], expected_hp[:, channel][keep]; rtol = 1.0e-11)
                    end
                end
            end
        end
    end
end
