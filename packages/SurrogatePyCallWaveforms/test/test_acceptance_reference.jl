using TestItems

@testsnippet Setup_acceptance_reference begin
    using Test
    using FFTW
    using GWSeries
    import AbstractWaveforms
    import SurrogateWaveforms
    import SurrogatePyCallWaveforms
    using AbstractWaveformsTests

    # The native surrogates against upstream gwsurrogate (through PyCall) as the reference generator.
    # Local only: needs the pinned Python environment and the model data, see
    # packages/SurrogatePyCallWaveforms/README.md and `just test-reference-surrogates`.
    const F_REF = 25.0
    # Rows inside the surrogates' calibrated range (mass ratio <= 4, stellar-mass BBH).
    const IDS = ("bbh_equal", "bbh_unequal", "bbh_aligned", "bbh_antialigned", "bbh_precessing")
    # Peak-relative times: the inspiral's last 0.1 s, merger and a short ringdown.
    const TIMES = range(-0.1, 0.01; step = 1 / 2048)
    const WINDOW = PlanckWindow(; nsamples = length(TIMES), sampling_interval = step(TIMES), start_width = 0.02, end_width = 0.005)
end

@testitem "SurrogateWaveforms against gwsurrogate" tags = [:integration, :reference, :python, :slow] setup = [Setup_acceptance_reference] begin
    for (model, ell_max) in ((:NRSur7dq4, 4), (:NRSur7dq4v2, 5))
        @testset "$model" begin
            surrogate = SurrogateWaveforms.load_surrogate(model)
            native = SurrogateWaveforms.PrecessingGenerator(surrogate; ell_max, f_ref = F_REF)
            reference = SurrogatePyCallWaveforms.load_generator(SurrogateWaveforms.data_path(model); model, ell_max, f_ref = F_REF)

            # Signed modes, on the model grid, a fixed time step and explicit times.
            cases = applicable_cases(native, reference; ids = IDS)
            labels = Tuple((l, m) for l in 2:ell_max for m in (-l):l)
            for (id, case) in cases
                params = surrogate_parameters(case)
                @testset "modes: $id" begin
                    for grid in ((;), (; dt = 1 / 2048), (; times = TIMES))
                        expected = AbstractWaveforms.waveform_modes_td(reference, params; grid...)
                        actual = AbstractWaveforms.waveform_modes_td(native, params; grid...)
                        @test actual.t ≈ expected.t
                        @test actual.mode_indices == expected.mode_indices == labels
                        # Upstream's C++ spline mis-evaluates isolated points in the last coorbital
                        # interval; shape and finiteness are still checked for every sample.
                        tscale = params.total_mass * SurrogateWaveforms.MSUN_IN_SECONDS
                        keep = expected.t .<= surrogate.coorbital.t[end - 1] * tscale
                        @test any(keep)
                        AbstractWaveformsTests.check_pointwise(actual.modes[keep, :], expected.modes[keep, :]; rtol = 1.0e-11)
                    end
                end
            end

            # Windowed polarizations compared through the physical rfft.
            test_reference(
                native, reference, cases; grid = TIMES, parameters = surrogate_parameters, target_domain = :td,
                prepare = h -> apply_window(WINDOW, h), max_mismatch = 1.0e-14, rtol = 1.0e-9,
            )
        end
    end
end
