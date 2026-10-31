# Shared optimal-SNR setup for tagged package acceptance tests. Needs GWNoise and
# GWMatchedFiltering in the including environment.
using GWNoise, GWMatchedFiltering

const PSD_NAMES = (:ET_000A_18_ETD_hlf, :cosmic_explorer_strain)
# A fixed response with both polarizations observable; geometry is not differentiated.
project_strain(h) = 0.7 .* h[:, 1] .+ 0.3 .* h[:, 2]

# One-sided PSDs in strain²/Hz read through GWNoise with log-log interpolation. The analysis band is
# 20–512 Hz intersected with table and grid coverage; only active bins are interpolated, and bins
# outside the band, including DC, are zero.
#
# Returns the `(nf, 1)` whitening factor; the caller passes `step(f)` alongside it.
function snr_whitening(f, name)
    curve = TabulatedPSD(psd_path(name), ' '; interpolation = :loglog)
    low = max(20.0, first(curve.frequencies))
    high = min(512.0, last(curve.frequencies), last(f))
    return whitening_factor(reshape(band_psd(curve, f; fmin = low, fmax = high), :, 1))
end
