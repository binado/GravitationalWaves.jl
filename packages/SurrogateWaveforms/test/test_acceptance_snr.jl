using TestItems

@testsnippet Setup_acceptance_snr begin
    using GWSeries
    using Test
    import AbstractWaveforms
    using AbstractWaveforms: mode_sum, sylm_m2, HigherModesWaveformGenerator, waveform_modes_fd, waveform_modes_td, hphc_fd, hphc_td
    using DifferentiationInterface
    using FiniteDifferences: central_fdm, grad
    using ForwardDiff: ForwardDiff
    using SurrogateWaveforms
    using Mooncake: Mooncake

    using AbstractWaveforms: hphc_td
    using ADTypes: AutoForwardDiff, AutoMooncake
    using AbstractWaveformsTests
    using FFTW, ChainRulesCore
    include(joinpath(pkgdir(AbstractWaveformsTests), "test", "support", "snr.jl"))

    const SNR_IDS = ("bbh_equal", "bbh_unequal", "bbh_aligned", "bbh_antialigned", "bbh_precessing")
    const SNR_BACKENDS = (AutoForwardDiff(), AutoMooncake(; config = nothing))

    # Fixed configuration is global, not captured: Mooncake would otherwise try to
    # copy the TwicePrecision range metadata.
    # Fixed peak-relative physical grid, including merger and a short ringdown.
    const SNR_TIMES = range(-0.2; step = 1 / 2048, length = 432)
    const SNR_WINDOW = PlanckWindow(; nsamples = length(SNR_TIMES), sampling_interval = step(SNR_TIMES), start_width = 0.02, end_width = 0.005)
    const SNR_PADDING = zeros(80, 2)
    const SNR_META = (; epoch = first(SNR_TIMES), sampling_interval = step(SNR_TIMES), channels = (:plus, :cross))
end

@testitem "TD optimal SNR" tags = [:integration, :ad, :slow] setup = [Setup_acceptance_snr] begin
    t = SNR_TIMES
    n = length(t) + size(SNR_PADDING, 1)
    f = LinRange(0.0, (n ÷ 2) / (n * step(t)), n ÷ 2 + 1)
    for model in (:NRSur7dq4, :NRSur7dq4v2)
        sur = SurrogateWaveforms.load_surrogate(model)
        # f_ref=0 fixes the reference to the first dynamics node, avoiding a
        # parameter-dependent frequency-to-reference-time search.
        generator = SurrogateWaveforms.PrecessingGenerator(sur; f_ref = 0.0)
        for (id, p) in applicable_cases(generator; ids = SNR_IDS)
            @testset "$model / $id" begin
                x, scales = scaled_coordinates(p; precessing = true)
                adapter = surrogate_adapter(scales)
                function strain(z)
                    h = hphc_td(SNR_TIMES, adapter(z), generator)
                    samples = vcat(apply_window(SNR_WINDOW, h), SNR_PADDING)
                    spectrum = to_frequencydomain(TimeSeries(samples; SNR_META...))
                    return project_strain(parent(spectrum))
                end
                for name in PSD_NAMES
                    @testset "$name" begin
                        W = snr_whitening(f, name)
                        test_snr_gradients(strain, W, step(f), x; backends = SNR_BACKENDS, rtol = 1.0e-6, finite_difference = id == "bbh_precessing")
                    end
                end
                id == "bbh_precessing" && test_phase_gradients(strain, x; backends = SNR_BACKENDS, rtol = 1.0e-6)
            end
        end
    end
end
