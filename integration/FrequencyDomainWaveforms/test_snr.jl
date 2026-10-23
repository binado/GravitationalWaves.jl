using ADTypes: AutoForwardDiff, AutoMooncake
using AbstractWaveformsTests
using Mooncake: Mooncake
include("../common/snr.jl")

const SNR_GRID = LinRange(0.0, 1024.0, 2049)
const SNR_IDS = ("bbh_equal", "bbh_unequal", "bbh_aligned", "bbh_antialigned")
const SNR_BACKENDS = (AutoForwardDiff(), AutoMooncake(; config = nothing))

@testset "FD optimal SNR" begin
    f = SNR_GRID
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        generator = G(; minimum_frequency = 20.0, maximum_frequency = 1024.0)
        for (id, p) in applicable_cases(generator; ids = SNR_IDS)
            @testset "$G / $id" begin
                x, scales = scaled_coordinates(p)
                adapter = fd_adapter(p, scales)
                validate_parameters(f, adapter(x), generator)
                strain(z) = project_strain(hphc_fd(SNR_GRID, adapter(z), generator))
                for name in PSD_NAMES
                    @testset "$name" begin
                        metric = snr_metric(f, name)
                        test_snr_gradients(strain, metric, x; backends = SNR_BACKENDS, rtol = 1.0e-8, finite_difference = id == "bbh_aligned")
                    end
                end
                id == "bbh_aligned" && test_phase_gradients(strain, x; backends = SNR_BACKENDS, rtol = 1.0e-6)
            end
        end
    end
end

@testset "FD waveform Jacobian" begin
    for G in (TaylorF2, IMRPhenomD, IMRPhenomXAS)
        generator = G(; minimum_frequency = 20.0, maximum_frequency = 1024.0)
        @testset "$G" begin
            p = default_cases()["bbh_aligned"]
            x, scales = scaled_coordinates(p)
            # Forward mode only: the output has thousands of components.
            test_waveform_jacobian(
                generator, x, fd_adapter(p, scales); grid = 0.0:4.0:1024.0, backends = (AutoForwardDiff(),), rtol = 1.0e-8,
                fd_rtol = 1.0e-5,
            )
        end
    end
end
