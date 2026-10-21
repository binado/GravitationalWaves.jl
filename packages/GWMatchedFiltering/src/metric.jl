function analysis_bins(f, minimum_frequency, maximum_frequency, nsamples)
    isfinite(minimum_frequency) && isfinite(maximum_frequency) &&
        first(f) <= minimum_frequency <= maximum_frequency <= last(f) ||
        throw(ArgumentError("analysis band must lie within the supplied grid"))
    lower = max(searchsortedfirst(f, minimum_frequency), searchsortedlast(f, 0) + 1)
    upper = searchsortedlast(f, maximum_frequency)
    if !isnothing(nsamples)
        nsamples isa Integer && nsamples >= 2 || throw(ArgumentError("nsamples must be an integer >= 2"))
        first(f) == 0 && length(f) == nsamples ÷ 2 + 1 ||
            throw(ArgumentError("nsamples requires a complete real-FFT frequency grid starting at zero"))
        iseven(nsamples) && (upper = min(upper, length(f) - 1))
    end
    lower <= upper || throw(ArgumentError("analysis band contains no active interior bins"))
    return lower:upper
end

"""
    FrequencyDomainMetric(series::FrequencySeries, spectral_density;
        minimum_frequency=first(frequencies(series)),
        maximum_frequency=last(frequencies(series)), nsamples=nothing)

Prepare a joint metric with first-argument conjugation and normalization `4Δf` stored through the metric frequency grid.
The domain frequency grid must be a uniform nonnegative `AbstractRange`. Supply PSDs with shape
`(frequency, channel)` (or a vector for one channel), or Hermitian spectral
matrices with shape `(frequency, channel, channel)`. Define cross spectra by
`S[k,i,j] = 2Δf * E[n[k,i] * conj(n[k,j])]` for interior bins.

Band bounds are inclusive: first bin >= minimum, last bin <= maximum. DC is
always excluded. Supply `nsamples` for a complete real-FFT grid to also exclude
even-length Nyquist; an odd-length final bin remains active. Without `nsamples`,
the final positive-frequency bin is ordinary (Nyquist cannot be inferred).

Only active densities are inspected. Active PSDs must be finite and positive;
active spectral matrices must be finite, exactly Hermitian and positive definite.
Zeros outside the band are allowed. No interpolation, jitter, pseudoinverse or
automatic channel reordering is performed. Setup copies the prepared weights.
"""
function FrequencyDomainMetric(
        spectrum::AbstractArray{<:Number}; frequencies::AbstractRange, channels,
        minimum_frequency = first(frequencies), maximum_frequency = last(frequencies), nsamples = nothing,
    )
    check_grid(frequencies)
    f = range(GWSeries.metadata_value(first(frequencies)); step = GWSeries.metadata_value(step(frequencies)), length = length(frequencies))
    names = GWSeries.channel_order(channels)
    bins = analysis_bins(f, minimum_frequency, maximum_frequency, nsamples)
    prepared = prepare_metric(spectrum, length(f), length(names), bins)
    return FrequencyDomainMetric(f, names, prepared, bins)
end
function FrequencyDomainMetric(series::FrequencySeries, spectrum::AbstractArray{<:Number}; kwargs...)
    metadata(series).spectrum == :onesided || throw(ArgumentError("one-sided metrics require a one-sided spectrum"))
    return FrequencyDomainMetric(spectrum; frequencies = frequencies(series), channels = channels(series), kwargs...)
end
frequency_spacing(metric::FrequencyDomainMetric) = step(metric.frequencies)

prepare_metric(spectrum::AbstractVector{<:Real}, nf, nc, bins) =
    prepare_metric(reshape(spectrum, :, 1), nf, nc, bins)

function prepare_metric(spectrum::AbstractMatrix{<:Real}, nf, nc, bins)
    size(spectrum) == (nf, nc) || throw(DimensionMismatch("PSDs must have shape (frequency, channel)"))
    selected = spectrum[bins, :]
    all(x -> isfinite(x) && x > 0, selected) ||
        throw(ArgumentError("active PSDs must be finite and positive"))
    return DiagonalMetric(inv.(selected))
end

function prepare_metric(spectrum::AbstractArray{<:Number, 3}, nf, nc, bins)
    size(spectrum) == (nf, nc, nc) ||
        throw(DimensionMismatch("spectral matrices must have shape (frequency, channel, channel)"))
    # Host setup: validate before Hermitian(), which otherwise trusts one triangle.
    T = typeof(float(zero(eltype(spectrum))))
    factors = Array{T}(undef, nc, nc, length(bins))
    for (k, source) in enumerate(bins)
        block = Matrix{T}(spectrum[source, :, :])
        all(isfinite, block) && ishermitian(block) ||
            throw(ArgumentError("active spectral matrices must be finite and Hermitian"))
        factor = cholesky(Hermitian(block))
        factors[:, :, k] = factor.L \ Matrix{T}(I, nc, nc)
    end
    return BlockDiagonalMetric(factors)
end

"""Host-only stationary metric compatibility; the segment epoch is irrelevant."""
function check_compatible(d::FrequencySeries, metric::FrequencyDomainMetric)
    metadata(d).spectrum == :onesided || throw(ArgumentError("one-sided metrics require a one-sided spectrum"))
    check_same_frequency_grid(frequencies(d), metric.frequencies)
    channels(d) == metric.channels || throw(ArgumentError("metric channel orders differ"))
    return nothing
end
