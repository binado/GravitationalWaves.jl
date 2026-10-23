# Matched-filtering operations on plain arrays and a whitening factor `W`.
#
# Conventions (defined in GWNoise): data are `(nf, nc)` complex one-sided
# frequency-domain arrays (a single channel may be an `(nf,)` vector); `W` is
# `L⁻¹` with `L L† = S` per bin, either `(nf, nc)` real (diagonal, `1/√S`, or an
# `(nf,)` vector for one channel) or `(nc, nc, nf)` lower triangular (blocks).
# Bins where `W` is zero are inactive. Kernels check no shapes or values.

# Diagonal weights for a single-channel vector accept `(nf,)` or `(nf, 1)` storage.
diagonal_weight(W::AbstractVector, k) = @inbounds W[k]
diagonal_weight(W::AbstractMatrix, k) = @inbounds W[k, 1]

# Reactant's traced arrays also subtype DenseArray/StridedArray. Restrict the
# native specializations to host storage (and its views), so device arrays keep
# frequency-wise broadcast/reduction rather than unrolling scalar loops.
const HostMatrix = Union{Matrix, SubArray{<:Any, 2, <:Matrix}}
# A vector is one channel. Keep it one-dimensional: reshaping constant input
# buffers inside a differentiated call introduces mixed-activity array headers.
const HostVector = Union{Vector, SubArray{<:Any, 1, <:Union{Vector, Matrix}}}

# Σₖ conj(a) W² b over independent channels.
whitened_contraction(a, b, W::AbstractMatrix) = sum(conj.(a) .* abs2.(W) .* b)
function whitened_contraction(a::HostMatrix, b::HostMatrix, W::HostMatrix)
    result = zero(promote_type(eltype(a), eltype(b), eltype(W)))
    @inbounds for channel in axes(a, 2), k in axes(a, 1)
        w = W[k, channel]
        result += conj(a[k, channel]) * (w * w) * b[k, channel]
    end
    return result
end

function whitened_contraction(a::AbstractVector, b::AbstractVector, W::AbstractMatrix)
    return sum(conj.(a) .* abs2.(view(W, :, 1)) .* b)
end
function whitened_contraction(a::AbstractVector, b::AbstractVector, W::AbstractVector)
    return sum(conj.(a) .* abs2.(W) .* b)
end
function whitened_contraction(a::HostVector, b::HostVector, W::Union{HostVector, HostMatrix})
    result = zero(promote_type(eltype(a), eltype(b), eltype(W)))
    @inbounds for k in eachindex(a, b)
        w = diagonal_weight(W, k)
        result += conj(a[k]) * (w * w) * b[k]
    end
    return result
end

# Σₖ (W aₖ)† (W bₖ) for `(nc, nc, nf)` blocks.
const HostArray3 = Union{Array{<:Any, 3}, SubArray{<:Any, 3, <:Array}}

function whitened_contraction(a::HostMatrix, b::HostMatrix, W::HostArray3)
    result = zero(promote_type(eltype(a), eltype(b), eltype(W)))
    @inbounds for k in axes(a, 1), i in axes(W, 1)
        wa = zero(result)
        wb = zero(result)
        for j in axes(W, 2)
            wa += W[i, j, k] * a[k, j]
            wb += W[i, j, k] * b[k, j]
        end
        result += conj(wa) * wb
    end
    return result
end

# A channel-sized host loop, with array operations over frequencies. This
# avoids mutation and supports compilation/differentiation with a fixed `W`.
function whitened_contraction(a::AbstractMatrix, b::AbstractMatrix, W::AbstractArray{<:Any, 3})
    result = nothing
    for channel in axes(W, 1)
        weights = permutedims(W[channel, :, :], (2, 1))
        white_a = sum(weights .* a; dims = 2)
        white_b = sum(weights .* b; dims = 2)
        term = sum(conj.(white_a) .* white_b)
        result = isnothing(result) ? term : result + term
    end
    return result
end

# One channel stored in `(1, 1, nf)` blocks.
function whitened_contraction(a::AbstractVector, b::AbstractVector, W::AbstractArray{<:Any, 3})
    weights = view(W, 1, 1, :)
    return sum(conj.(weights .* a) .* (weights .* b))
end
function whitened_contraction(a::HostVector, b::HostVector, W::HostArray3)
    result = zero(promote_type(eltype(a), eltype(b), eltype(W)))
    @inbounds for k in eachindex(a, b)
        w = W[1, 1, k]
        result += conj(w * a[k]) * (w * b[k])
    end
    return result
end

"""
    complex_inner_product(a, b, W, df)

Compute `4df Σₖ (W aₖ)† (W bₖ)` for data in `(frequency, channel)` layout (a vector
is one channel) and a whitening factor `W`: `(frequency, channel)` real for
independent channels, `(channel, channel, frequency)` for correlated ones. The
first argument is conjugated. Bins where `W` is zero drop out. Arrays are
prevalidated; no shape or value checks run.
"""
complex_inner_product(a::AbstractMatrix, b::AbstractMatrix, W::AbstractArray, df) =
    (4 * df) * whitened_contraction(a, b, W)
complex_inner_product(a::AbstractVector, b::AbstractVector, W::AbstractArray, df) =
    (4 * df) * whitened_contraction(a, b, W)
complex_inner_product(a::AbstractVector, b::AbstractMatrix, W::AbstractArray, df) =
    complex_inner_product(a, view(b, :, 1), W, df)
complex_inner_product(a::AbstractMatrix, b::AbstractVector, W::AbstractArray, df) =
    complex_inner_product(view(a, :, 1), b, W, df)

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
`W[:, :, k] * d[k, :]` for blocks. `out` must not alias `d` for blocks.
"""
whiten(d::AbstractVecOrMat, W::AbstractVecOrMat) = diagonal_weights(W, d) .* d
diagonal_weights(W::AbstractVector, ::AbstractVector) = W
diagonal_weights(W::AbstractMatrix, ::AbstractMatrix) = W
diagonal_weights(W::AbstractMatrix, ::AbstractVector) = view(W, :, 1)
function whiten(d::AbstractMatrix, W::AbstractArray{<:Any, 3})
    return whiten!(similar(d, promote_type(eltype(d), eltype(W))), d, W)
end
whiten!(out::AbstractVecOrMat, d::AbstractVecOrMat, W::AbstractVecOrMat) = (out .= whiten(d, W))
function whiten!(out::AbstractMatrix, d::AbstractMatrix, W::AbstractArray{<:Any, 3})
    @inbounds for k in axes(d, 1), i in axes(W, 1)
        acc = zero(eltype(out))
        for j in axes(W, 2)
            acc += W[i, j, k] * d[k, j]
        end
        out[k, i] = acc
    end
    return out
end
