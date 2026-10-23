"""
    flat_metric(f; channels = (:plus, :cross), minimum_frequency = first(f), maximum_frequency = last(f), nsamples = nothing)

Unit-PSD `FrequencyDomainMetric` on the uniform grid `f`, so inner products are plain sums.
"""
function flat_metric(
        f::AbstractRange; channels = (:plus, :cross),
        minimum_frequency = first(f), maximum_frequency = last(f), nsamples = nothing,
    )
    return FrequencyDomainMetric(
        ones(length(f), length(channels)); frequencies = f, channels,
        minimum_frequency, maximum_frequency, nsamples,
    )
end

"""
    mismatch(a, b, metric; maximize = :time_phase, max_shift = 0.05)

`1 - overlap` of two `(frequency, channel)` polarization arrays (or `FrequencySeries`) on the
grid of `metric`, with the noise-weighted inner product of GWMatchedFiltering.

`maximize` selects what is optimized over a common shift of all channels:
- `:none`: the real part of the overlap; any time or phase offset counts as mismatch.
- `:phase`: the modulus of the complex overlap, i.e. a constant phase.
- `:time_phase`: the modulus, also maximized over a time shift in `[-max_shift, max_shift]` seconds
  (grid scan at 1/8 of the shortest period, then golden-section refinement).

The result is clamped at zero, since rounding can push the overlap slightly above one.
"""
function mismatch(a::AbstractMatrix, b::AbstractMatrix, metric::FrequencyDomainMetric; maximize::Symbol = :time_phase, max_shift = 0.05)
    maximize in (:none, :phase, :time_phase) || throw(ArgumentError("maximize must be :none, :phase or :time_phase"))
    size(a) == size(b) || throw(DimensionMismatch("waveforms differ in shape: $(size(a)) vs $(size(b))"))
    norm_a = optimal_snr_squared(a, metric)
    norm_b = optimal_snr_squared(b, metric)
    (isfinite(norm_a) && norm_a > 0 && isfinite(norm_b) && norm_b > 0) ||
        throw(DomainError((norm_a, norm_b), "mismatch needs finite nonzero waveform norms in the metric band"))
    scale = sqrt(norm_a * norm_b)
    if maximize === :none
        overlap = real(complex_inner_product(a, b, metric)) / scale
    elseif maximize === :phase
        overlap = abs(complex_inner_product(a, b, metric)) / scale
    else
        overlap = maximum_shifted_overlap(a, b, metric, scale, max_shift)
    end
    return max(zero(overlap), 1 - overlap)
end
mismatch(a::FrequencySeries, b::FrequencySeries, metric::FrequencyDomainMetric; kwargs...) =
    mismatch(parent(a), parent(b), metric; kwargs...)

# The metric's grid is a stored field; there is no accessor for it.
function maximum_shifted_overlap(a, b, metric, scale, max_shift)
    f = collect(float.(metric.frequencies))
    at(τ) = abs(complex_inner_product(a, b .* cis.((-2π * τ) .* f), metric)) / scale
    step_τ = 1 / (8 * maximum(f))
    taus = (-max_shift):step_τ:max_shift
    values = map(at, taus)
    best = argmax(values)
    lo, hi = taus[max(best - 1, 1)], taus[min(best + 1, length(taus))]
    # Golden-section search; the overlap is smooth and unimodal within one grid step of the peak.
    ϕ = (sqrt(5) - 1) / 2
    x1, x2 = hi - ϕ * (hi - lo), lo + ϕ * (hi - lo)
    f1, f2 = at(x1), at(x2)
    for _ in 1:60
        if f1 < f2
            lo, x1, f1 = x1, x2, f2
            x2 = lo + ϕ * (hi - lo)
            f2 = at(x2)
        else
            hi, x2, f2 = x2, x1, f1
            x1 = hi - ϕ * (hi - lo)
            f1 = at(x1)
        end
    end
    return max(values[best], f1, f2)
end
