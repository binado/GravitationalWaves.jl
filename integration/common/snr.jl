# Shared optimal-SNR setup for the heavy gradient suites under integration/. Needs GWNoise and
# GWMatchedFiltering in the including environment.
using GWNoise, GWMatchedFiltering

const PSD_NAMES = (:ET_000A_18_ETD_hlf, :cosmic_explorer_strain)
# A fixed response with both polarizations observable; geometry is not differentiated.
project_strain(h) = 0.7 .* h[:, 1] .+ 0.3 .* h[:, 2]

# One-sided PSDs in strain²/Hz read through GWNoise with log-log interpolation. The analysis band is
# 20–512 Hz intersected with table and grid coverage; only active bins are interpolated.
function snr_metric(f, name; nsamples = nothing)
    curve = read_psd(psd_path(name); interpolation = :loglog)
    low = max(20.0, first(curve.frequencies))
    high = min(512.0, last(curve.frequencies), last(f))
    spectrum = [low <= x <= high ? psd(x, curve) : 0.0 for x in f]
    return FrequencyDomainMetric(
        spectrum; frequencies = f, channels = (:strain,),
        minimum_frequency = low, maximum_frequency = high, nsamples
    )
end
