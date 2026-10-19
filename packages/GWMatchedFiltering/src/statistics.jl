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
    complex_inner_product(domain, a, b, metric)

Compute `4df Σ a[k,:]' inv(S[k,:,:]) b[k,:]`. Arrays use (frequency, channel)
layout; vectors denote one channel. Validate data/domain/metric compatibility on
the host before evaluation. No metadata or sample-value checks occur here.
"""
function complex_inner_product(d::FrequencyDomain, a::AbstractMatrix, b::AbstractMatrix, metric::FrequencyDomainMetric)
    return (4 * frequency_spacing(d)) * contraction(active_values(a, metric.active_bins), active_values(b, metric.active_bins), metric.metric)
end
"""Real noise-weighted inner product of prevalidated arrays."""
noise_weighted_inner_product(d::FrequencyDomain, a, b, metric::FrequencyDomainMetric) = real(complex_inner_product(d, a, b, metric))
"""Squared optimal SNR for prevalidated samples and fixed sampling/metric metadata."""
optimal_snr_squared(d::FrequencyDomain, h, metric::FrequencyDomainMetric) = noise_weighted_inner_product(d, h, h, metric)
"""Optimal SNR for prevalidated arrays."""
optimal_snr(d::FrequencyDomain, h, metric::FrequencyDomainMetric) = sqrt(optimal_snr_squared(d, h, metric))
"""Complex scalar matched-filter SNR. Requires a nonzero template; run check_template on the host."""
matched_filter_snr(d::FrequencyDomain, template, data, metric::FrequencyDomainMetric) =
    complex_inner_product(d, template, data, metric) / sqrt(optimal_snr_squared(d, template, metric))
"""Squared residual norm for prevalidated arrays."""
function residual_squared(d::FrequencyDomain, data::AbstractMatrix, template::AbstractMatrix, metric::FrequencyDomainMetric)
    residual = data .- template
    return optimal_snr_squared(d, residual, metric)
end
"""Reusable fixed-PSD contractions, without metadata validation."""
likelihood_statistics(d::FrequencyDomain, template, data, metric::FrequencyDomainMetric) =
    (; h_inner_d = complex_inner_product(d, template, data, metric), h_inner_h = optimal_snr_squared(d, template, metric))
"""Fixed-PSD Gaussian log likelihood ratio for prevalidated arrays."""
function log_likelihood_ratio(d::FrequencyDomain, template, data, metric::FrequencyDomainMetric)
    stats = likelihood_statistics(d, template, data, metric)
    return real(stats.h_inner_d) - stats.h_inner_h / 2
end
"""Host-only check of domain compatibility, data shape, and finite positive template norm."""
function check_template(d::FrequencyDomain, template, metric::FrequencyDomainMetric)
    check_compatible(d, metric)
    check_data(d, template)
    norm_squared = optimal_snr_squared(d, template, metric)
    isfinite(norm_squared) && norm_squared > 0 || throw(DomainError(norm_squared, "matched-filter SNR needs a finite nonzero template"))
    return nothing
end

# A vector is one channel. Keep it one-dimensional: reshaping constant input
# buffers inside a differentiated call introduces mixed-activity array headers.
const HostVector = Union{Vector, SubArray{<:Any, 1, <:Union{Vector, Matrix}}}
active_values(a::AbstractVector, bins) = a[bins]
# Index original arrays directly; Base.view(::Vector, ...) itself invokes an
# identity reshape that Enzyme cannot share between active and constant buffers.
function contraction(a::HostVector, b::HostVector, metric::DiagonalMetric, bins)
    result = zero(promote_type(eltype(a), eltype(b), eltype(metric.weights)))
    @inbounds for (j, k) in enumerate(bins)
        result += conj(a[k]) * metric.weights[j, 1] * b[k]
    end
    return result
end
function contraction(a::HostVector, b::HostVector, metric::BlockDiagonalMetric{1}, bins)
    result = zero(promote_type(eltype(a), eltype(b), eltype(metric.whiteners)))
    @inbounds for (j, k) in enumerate(bins)
        w = metric.whiteners[1, 1, j]
        result += conj(w * a[k]) * (w * b[k])
    end
    return result
end
function complex_inner_product(d::FrequencyDomain, a::HostVector, b::HostVector, m::FrequencyDomainMetric)
    return (4 * frequency_spacing(d)) * contraction(a, b, m.metric, m.active_bins)
end
contraction(a::AbstractVector, b::AbstractVector, metric::DiagonalMetric) =
    sum(conj.(a) .* view(metric.weights, :, 1) .* b)
function contraction(a::AbstractVector, b::AbstractVector, metric::BlockDiagonalMetric{1})
    weights = view(metric.whiteners, 1, 1, :)
    return sum(conj.(weights .* a) .* (weights .* b))
end
function complex_inner_product(d::FrequencyDomain, a::AbstractVector, b::AbstractVector, m::FrequencyDomainMetric)
    return (4 * frequency_spacing(d)) * contraction(active_values(a, m.active_bins), active_values(b, m.active_bins), m.metric)
end
complex_inner_product(d::FrequencyDomain, a::AbstractVector, b::AbstractMatrix, m::FrequencyDomainMetric) =
    complex_inner_product(d, a, view(b, :, 1), m)
complex_inner_product(d::FrequencyDomain, a::AbstractMatrix, b::AbstractVector, m::FrequencyDomainMetric) =
    complex_inner_product(d, view(a, :, 1), b, m)
residual_squared(d::FrequencyDomain, a::AbstractVector, b::AbstractVector, m::FrequencyDomainMetric) =
    optimal_snr_squared(d, a .- b, m)
residual_squared(d::FrequencyDomain, a::AbstractVector, b::AbstractMatrix, m::FrequencyDomainMetric) =
    residual_squared(d, a, view(b, :, 1), m)
residual_squared(d::FrequencyDomain, a::AbstractMatrix, b::AbstractVector, m::FrequencyDomainMetric) =
    residual_squared(d, view(a, :, 1), b, m)
