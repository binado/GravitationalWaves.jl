using AbstractWaveforms: hphc_td
using FFTW, ChainRulesCore
include("../../test/waveforms/cases.jl")
include("../../test/waveforms/ad_checks.jl")

# Fixed configuration is global, not captured: Mooncake would otherwise try to
# copy the TwicePrecision range metadata.
# Fixed peak-relative physical grid, including merger and a short ringdown.
const SNR_TIMES = range(-0.2; step = 1 / 2048, length = 432)
const SNR_WINDOW = PlanckWindow(; nsamples = length(SNR_TIMES), sampling_interval = step(SNR_TIMES), start_width = 0.02, end_width = 0.005)
const SNR_PADDING = zeros(80, 2)
const SNR_META = (; epoch = first(SNR_TIMES), sampling_interval = step(SNR_TIMES), channels = (:plus, :cross))

@testset "TD optimal SNR" begin
    t = SNR_TIMES
    n = length(t) + size(SNR_PADDING, 1)
    f = LinRange(0.0, (n ÷ 2) / (n * step(t)), n ÷ 2 + 1)
    for model in (:NRSur7dq4, :NRSur7dq4v2)
        sur = SurrogateWaveforms.load_surrogate(model)
        # f_ref=0 fixes the reference to the first dynamics node, avoiding a
        # parameter-dependent frequency-to-reference-time search.
        generator = SurrogateWaveforms.PrecessingGenerator(sur; f_ref = 0.0)
        for id in WaveformCases.SURROGATE_IDS
            @testset "$model / $id" begin
                p = WaveformCases.CASES[id]
                x, scales = WaveformCases.scaled_coordinates(p; precessing = true)
                adapter = WaveformCases.surrogate_adapter(scales)
                function strain(z)
                    h = hphc_td(SNR_TIMES, adapter(z), generator)
                    samples = vcat(apply_window(SNR_WINDOW, h), SNR_PADDING)
                    spectrum = to_frequencydomain(TimeSeries(samples; SNR_META...))
                    return project_strain(parent(spectrum))
                end
                for name in PSD_NAMES
                    @testset "$name" begin
                        metric = snr_metric(f, name; nsamples = n)
                        loss(z) = sqrt(optimal_snr_squared(strain(z), metric))
                        check_snr_ad(loss, x; rtol = 1.0e-6, finite_difference = id == "bbh_precessing")
                    end
                end
                id == "bbh_precessing" && check_phase_ad(strain, x)
            end
        end
    end
end
