# Matched-filtering operations on plain arrays and a whitening factor `W`.
#
# Conventions (defined in GWNoise): data are `(nf, nc)` complex one-sided
# frequency-domain arrays (a single channel may be an `(nf,)` vector); `W` is
# `L⁻¹` with `L L† = S` per bin, either `(nf, nc)` real (diagonal, `1/√S`, or an
# `(nf,)` vector for one channel) or `(nc, nc, nf)` lower triangular (blocks).
# Bins where `W` is zero are inactive. Kernels check no shapes or values.

# How a whitening array `W` is read. Extend `filter_structure` to opt in a new storage form.
abstract type FilterStructure end
"""Independent channels: one real weight per bin and channel, `W = 1/√S`."""
struct DiagonalFilter <: FilterStructure end
"""Correlated channels: a lower-triangular `nc × nc` factor per bin, bins independent."""
struct BlockFilter <: FilterStructure end

filter_structure(::AbstractVector{<:Real}) = DiagonalFilter()
filter_structure(::AbstractMatrix{<:Real}) = DiagonalFilter()
filter_structure(::AbstractArray{<:Number, 3}) = BlockFilter()

# Σ (W aₖ)† (W bₖ). The diagonal case is a plain multi-array reduction (shapes
# match, or one side is a single-channel vector against an `(nf, 1)` array); blocks
# whiten both inputs with the Cholesky-derived factor and take an ordinary `dot`.
whitened_contraction(a, b, W) = whitened_contraction(filter_structure(W), a, b, W)
whitened_contraction(::DiagonalFilter, a, b, W) = mapreduce((x, w, y) -> conj(x) * (w * w) * y, +, a, W, b)
whitened_contraction(::BlockFilter, a, b, W) = dot(whiten(a, W), whiten(b, W))

"""
    complex_inner_product(a, b, W, df)

Compute `4df Σₖ (W aₖ)† (W bₖ)` for data in `(frequency, channel)` layout (a vector
is one channel) and a whitening factor `W`: `(frequency, channel)` real for
independent channels, `(channel, channel, frequency)` for correlated ones. The
first argument is conjugated. Bins where `W` is zero drop out. Arrays are
prevalidated; no shape or value checks run.
"""
complex_inner_product(a, b, W, df) = (4 * df) * whitened_contraction(a, b, W)

"""Real noise-weighted inner product, the real part of [`complex_inner_product`](@ref)."""
noise_weighted_inner_product(a, b, W, df) = real(complex_inner_product(a, b, W, df))
"""Squared optimal SNR `⟨h, h⟩`."""
optimal_snr_squared(h, W, df) = noise_weighted_inner_product(h, h, W, df)
"""Optimal SNR `√⟨h, h⟩`."""
optimal_snr(h, W, df) = sqrt(optimal_snr_squared(h, W, df))
"""Complex scalar matched-filter SNR. Requires a nonzero template; run `check_template` on the host."""
matched_filter_snr(template, data, W, df) =
    complex_inner_product(template, data, W, df) / sqrt(optimal_snr_squared(template, W, df))
"""Squared residual norm `⟨d - h, d - h⟩`."""
residual_squared(data, template, W, df) = optimal_snr_squared(data .- template, W, df)
"""Reusable fixed-PSD contractions `⟨h, d⟩` and `⟨h, h⟩`."""
likelihood_statistics(template, data, W, df) =
    (; h_inner_d = complex_inner_product(template, data, W, df), h_inner_h = optimal_snr_squared(template, W, df))
"""Fixed-PSD Gaussian log likelihood ratio `Re⟨h, d⟩ - ⟨h, h⟩/2`."""
function log_likelihood_ratio(template, data, W, df)
    stats = likelihood_statistics(template, data, W, df)
    return real(stats.h_inner_d) - stats.h_inner_h / 2
end

whitening_nbins(W::AbstractVecOrMat) = size(W, 1)
whitening_nbins(W::AbstractArray{<:Any, 3}) = size(W, 3)
whitening_nchannels(W::AbstractVector) = 1
whitening_nchannels(W::AbstractMatrix) = size(W, 2)
whitening_nchannels(W::AbstractArray{<:Any, 3}) = size(W, 1)

"""Host-only check of template shape against `W` and of a finite nonzero template norm."""
function check_template(template, W, df)
    Base.require_one_based_indexing(template)
    size(template, 1) == whitening_nbins(W) && size(template, 2) == whitening_nchannels(W) ||
        throw(DimensionMismatch("template shape differs from whitening factor"))
    norm_squared = optimal_snr_squared(template, W, df)
    isfinite(norm_squared) && norm_squared > 0 ||
        throw(DomainError(norm_squared, "matched-filter SNR needs a finite nonzero template"))
    return nothing
end

"""
    whiten(d, W)
    whiten!(out, d, W)

Apply the whitening factor per frequency bin: `W[k, c] * d[k, c]` for diagonal `W`,
`W[:, :, k] * d[k, :]` for blocks. Written as array expressions, so it composes with
automatic differentiation; `out` may not alias `d`.
"""
whiten(d, W) = whiten(filter_structure(W), d, W)
whiten!(out, d, W) = copyto!(out, whiten(d, W))
whiten(::DiagonalFilter, d, W) = reshape(W, size(d)) .* d
# `Wp[k, i, j] * d[k, 1, j]` summed over `j`.
function whiten(::BlockFilter, d::AbstractMatrix, W)
    return dropdims(sum(permutedims(W, (3, 1, 2)) .* reshape(d, size(d, 1), 1, :); dims = 3); dims = 3)
end
whiten(::BlockFilter, d::AbstractVector, W) = view(W, 1, 1, :) .* d
