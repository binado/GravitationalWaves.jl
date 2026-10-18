function check_array(a::AbstractMatrix, metric::FrequencyDomainMetric)
    size(a) == (length(metric.frequencies), length(metric.channels)) ||
        throw(DimensionMismatch("strain must have shape (frequency, channel) matching the metric"))
    return a
end

check_array(a::AbstractVector, metric::FrequencyDomainMetric) = check_array(reshape(a, :, 1), metric)

function check_array(a::FrequencyDomainStrain, metric::FrequencyDomainMetric)
    a.frequencies == metric.frequencies || throw(ArgumentError("strain and metric frequency grids differ"))
    a.channels == metric.channels || throw(ArgumentError("strain and metric channel orders differ"))
    return check_array(a.values, metric)
end

contraction(a, b, metric::DiagonalMetric) = sum(conj.(a) .* metric.weights .* b)

# Reactant's traced arrays also subtype DenseArray/StridedArray. Restrict the
# native specialization to host Matrix storage (and its views), so device arrays
# retain frequency-wise broadcast/reduction rather than unrolling scalar loops.
const HostMatrix = Union{Matrix, SubArray{<:Any, 2, <:Matrix}}

function contraction(a::HostMatrix, b::HostMatrix, metric::DiagonalMetric)
    weights = metric.weights
    result = zero(promote_type(eltype(a), eltype(b), eltype(weights)))
    @inbounds for channel in axes(a, 2), k in axes(a, 1)
        result += conj(a[k, channel]) * weights[k, channel] * b[k, channel]
    end
    return result
end

function contraction(a::HostMatrix, b::HostMatrix, metric::BlockDiagonalMetric{D}) where {D}
    result = zero(promote_type(eltype(a), eltype(b), eltype(metric.whiteners)))
    @inbounds for k in axes(a, 1)
        aa = SVector{D}(ntuple(i -> a[k, i], Val(D)))
        bb = SVector{D}(ntuple(i -> b[k, i], Val(D)))
        W = SMatrix{D, D}(ntuple(n -> metric.whiteners[mod1(n, D), cld(n, D), k], Val(D * D)))
        wa, wb = W * aa, W * bb
        result += sum(conj.(wa) .* wb)
    end
    return result
end

active_values(a, bins) = a[bins, :]
active_values(a::HostMatrix, bins) = view(a, bins, :)

# A channel-sized static host loop, with array operations over frequencies. This
# avoids mutation and supports compilation/differentiation with a fixed metric.
function contraction(a, b, metric::BlockDiagonalMetric)
    first_weights = permutedims(metric.whiteners[1, :, :], (2, 1))
    first_a = sum(first_weights .* a; dims = 2)
    first_b = sum(first_weights .* b; dims = 2)
    result = sum(conj.(first_a) .* first_b)
    for channel in 2:size(a, 2)
        weights = permutedims(metric.whiteners[channel, :, :], (2, 1))
        white_a = sum(weights .* a; dims = 2)
        white_b = sum(weights .* b; dims = 2)
        result += sum(conj.(white_a) .* white_b)
    end
    return result
end

"""
    complex_inner_product(a, b, metric)

Compute `4Δf Σ a[k,:]' * inv(S[k,:,:]) * b[k,:]` on active bins.
Accept strain containers or arrays; raw arrays implicitly use the metric's grid
and channel order. A vector denotes one channel. The first argument is conjugated.
"""
function complex_inner_product(a, b, metric::FrequencyDomainMetric)
    aa = check_array(a, metric)
    bb = check_array(b, metric)
    return metric.normalization * contraction(active_values(aa, metric.active_bins), active_values(bb, metric.active_bins), metric.metric)
end

"""The real part of [`complex_inner_product`](@ref), summed over all channels."""
noise_weighted_inner_product(a, b, metric::FrequencyDomainMetric) = real(complex_inner_product(a, b, metric))

"""Squared optimal network SNR `(h|h)` under the supplied joint noise metric."""
optimal_snr_squared(h, metric::FrequencyDomainMetric) = noise_weighted_inner_product(h, h, metric)

"""Optimal network SNR, `sqrt(optimal_snr_squared(h, metric))`."""
optimal_snr(h, metric::FrequencyDomainMetric) = sqrt(optimal_snr_squared(h, metric))

"""
Complex scalar matched-filter SNR `C(template,data)/sqrt(C(template,template))`
at the supplied alignment. No time/phase maximization is performed. A zero-norm
template raises `DomainError`. Set `check=false` inside compiled kernels after
establishing a nonzero template on the host; this removes the value-dependent
exception path. Shape and metadata checks remain enabled.
"""
function matched_filter_snr(template, data, metric::FrequencyDomainMetric; check::Bool = true)
    norm_squared = optimal_snr_squared(template, metric)
    if check
        norm_squared > 0 || throw(DomainError(norm_squared, "matched-filter SNR needs a nonzero template"))
    end
    return complex_inner_product(template, data, metric) / sqrt(norm_squared)
end

"""Squared residual norm `(data-template|data-template)` on the active channels/bins."""
function residual_squared(data, template, metric::FrequencyDomainMetric)
    residual = check_array(data, metric) .- check_array(template, metric)
    return optimal_snr_squared(residual, metric)
end

"""Return reusable `(h_inner_d, h_inner_h)` network contractions for a fixed PSD."""
function likelihood_statistics(template, data, metric::FrequencyDomainMetric)
    return (;
        h_inner_d = complex_inner_product(template, data, metric),
        h_inner_h = optimal_snr_squared(template, metric),
    )
end

"""Fixed-PSD Gaussian log likelihood ratio `Re C(h,d) - C(h,h)/2` against noise alone."""
function log_likelihood_ratio(template, data, metric::FrequencyDomainMetric)
    stats = likelihood_statistics(template, data, metric)
    return real(stats.h_inner_d) - stats.h_inner_h / 2
end
