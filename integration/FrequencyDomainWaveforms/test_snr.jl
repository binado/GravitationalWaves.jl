include("../../test/waveforms/cases.jl")
include("../../test/waveforms/ad_checks.jl")

# Fixed 0:0.5:1024 Hz grid, kept global so it is not captured by the objective.
const SNR_GRID = LinRange(0.0, 1024.0, 2049)

@testset "FD optimal SNR" begin
    f = SNR_GRID
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        generator = G(; minimum_frequency = 20.0, maximum_frequency = 1024.0)
        for id in WaveformCases.COMMON_IDS
            @testset "$G / $id" begin
                p = WaveformCases.CASES[id]
                x, scales = WaveformCases.scaled_coordinates(p)
                adapter = WaveformCases.fd_adapter(p, scales)
                validate_parameters(f, adapter(x), generator)
                strain(z) = project_strain(hphc_fd(SNR_GRID, adapter(z), generator))
                for name in PSD_NAMES
                    @testset "$name" begin
                        metric = snr_metric(f, name)
                        loss(z) = sqrt(optimal_snr_squared(strain(z), metric))
                        check_snr_ad(loss, x; rtol = 1.0e-8, finite_difference = id == "bbh_aligned")
                    end
                end
                id == "bbh_aligned" && check_phase_ad(strain, x)
            end
        end
    end
end
