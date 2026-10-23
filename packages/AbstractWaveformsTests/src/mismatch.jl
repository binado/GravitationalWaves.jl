"""
    flat_whitening(f; nchannels = 2, minimum_frequency = first(f), maximum_frequency = last(f), nsamples = nothing)

Unit-PSD whitening factor `W` on the uniform grid `f`, so inner products are plain sums: ones with
shape `(length(f), nchannels)` inside the inclusive band and zero outside. DC is always zero; with `nsamples`
(the original transform length) an even length also zeroes the Nyquist bin.
"""
function flat_whitening(
        f::AbstractRange; nchannels::Integer = 2,
        minimum_frequency = first(f), maximum_frequency = last(f), nsamples = nothing,
    )
    W = [(minimum_frequency <= x <= maximum_frequency && x > 0) ? 1.0 : 0.0 for x in f, _ in 1:nchannels]
    nsamples !== nothing && iseven(nsamples) && (W[end, :] .= 0)
    return W
end

"""
    mismatch(a, b, W, f; maximize = :time_phase, max_shift = 0.05)

`1 - overlap` of two `(frequency, channel)` polarization arrays (or `FrequencySeries`) on the
uniform grid `f`, with the noise-weighted inner product of GWMatchedFiltering for the whitening
factor `W` (see `GWNoise.whitening_factor`; `flat_whitening` gives a unit one). Inactive bins
(`W == 0`) drop out. For `FrequencySeries` arguments, `W` alone is given and `f` comes from `a`.

`maximize` selects what is optimized over a common shift of all channels:
- `:none`: the real part of the overlap; any time or phase offset counts as mismatch.
- `:phase`: the modulus of the complex overlap, i.e. a constant phase.
- `:time_phase`: the modulus, also maximized over a time shift in `[-max_shift, max_shift]` seconds
  (grid scan at 1/8 of the shortest period, then golden-section refinement).

The result is clamped at zero, since rounding can push the overlap slightly above one.
"""
function mismatch(a::AbstractMatrix, b::AbstractMatrix, W::AbstractArray, f::AbstractRange; maximize::Symbol = :time_phase, max_shift = 0.05)
    maximize in (:none, :phase, :time_phase) || throw(ArgumentError("maximize must be :none, :phase or :time_phase"))
    size(a) == size(b) || throw(DimensionMismatch("waveforms differ in shape: $(size(a)) vs $(size(b))"))
    df = step(f)
    norm_a = optimal_snr_squared(a, W, df)
    norm_b = optimal_snr_squared(b, W, df)
    (isfinite(norm_a) && norm_a > 0 && isfinite(norm_b) && norm_b > 0) ||
        throw(DomainError((norm_a, norm_b), "mismatch needs finite nonzero waveform norms in the whitened band"))
    scale = sqrt(norm_a * norm_b)
    if maximize === :none
        overlap = real(complex_inner_product(a, b, W, df)) / scale
    elseif maximize === :phase
        overlap = abs(complex_inner_product(a, b, W, df)) / scale
    else
        overlap = maximum_shifted_overlap(a, b, W, f, scale, max_shift)
    end
    return max(zero(overlap), 1 - overlap)
end
mismatch(a::FrequencySeries, b::FrequencySeries, W::AbstractArray; kwargs...) =
    mismatch(parent(a), parent(b), W, frequencies(a); kwargs...)

function maximum_shifted_overlap(a, b, W, grid, scale, max_shift)
    f = collect(float.(grid))
    df = step(grid)
    at(τ) = abs(complex_inner_product(a, b .* cis.((-2π * τ) .* f), W, df)) / scale
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
