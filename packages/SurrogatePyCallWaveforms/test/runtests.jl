using Test
using PyCall
import AbstractWaveforms
import SurrogateWaveforms
using AbstractWaveformsTests: test_modes_td_contract
using SurrogatePyCallWaveforms

# Needs a Python with gwsurrogate (PyCall built against it) and the NRSur7dq4 model data. Without
# them the tests are skipped locally and fail on CI.
function python_ready()
    return try
        pyimport("gwsurrogate.surrogate")
        true
    catch err
        err isa InterruptException && rethrow()
        false
    end
end

@testset "SurrogatePyCallWaveforms" begin
    if !python_ready()
        msg = "PyCall cannot import gwsurrogate (see the package README); skipping SurrogatePyCallWaveforms tests"
        get(ENV, "CI", "") == "true" ? error(msg) : @warn msg
    else
        data_path(model) = SurrogateWaveforms.data_path(model)
        params = (
            total_mass = 65.0, luminosity_distance = 500.0, mass_ratio = 1.24, spin_1 = (0.3, 0.2, 0.4),
            spin_2 = (-0.1, 0.5, -0.3), inclination = 0.4, phi_ref = 1.3,
        )
        times = range(-0.1, 0.01; step = 1 / 2048)
        g = load_generator(data_path(:NRSur7dq4); model = :NRSur7dq4, ell_max = 4, f_ref = 25.0)

        @testset "accessors" begin
            @test AbstractWaveforms.approximant(g) == "NRSur7dq4"
            @test AbstractWaveforms.reference_frequency(g) == 25.0
            @test AbstractWaveforms.spin_support(g) === :precessing
            @test !AbstractWaveforms.is_differentiable(g)
            @test AbstractWaveforms.supports_mode_array(g)
        end

        test_modes_td_contract(g; params, times)

        @testset "units" begin
            base = waveform_modes_td(g, params; times)
            # Time scales with the total mass; the strain amplitude with mass over distance.
            # The reference frequency is in Hz, so halve it to keep the same dimensionless initial condition.
            heavier = merge(params, (; total_mass = 2params.total_mass))
            g_half = load_generator(data_path(:NRSur7dq4); model = :NRSur7dq4, ell_max = 4, f_ref = 12.5)
            scaled = waveform_modes_td(g_half, heavier; times = 2 .* times)
            @test scaled.t ≈ 2 .* base.t
            @test scaled.modes ≈ 2 .* base.modes rtol = 1.0e-10
            farther = waveform_modes_td(g, merge(params, (; luminosity_distance = 2params.luminosity_distance)); times)
            @test farther.modes ≈ base.modes ./ 2 rtol = 1.0e-12
            @test base.mode_indices == Tuple((l, m) for l in 2:4 for m in (-l):l)
        end

        @testset "argument validation" begin
            @test_throws ArgumentError waveform_modes_td(g, params; dt = 1 / 2048, times)
            @test_throws ArgumentError load_generator(data_path(:NRSur7dq4); model = :nonsense, f_ref = 25.0)
            @test_throws ArgumentError load_generator(data_path(:NRSur7dq4); ell_max = 9, f_ref = 25.0)
            @test_throws ArgumentError load_generator(data_path(:NRSur7dq4); f_ref = 10.0, f_low = 20.0)
        end
    end
end
