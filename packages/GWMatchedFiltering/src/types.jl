"""Interface for prepared noise-weighted contractions."""
abstract type AbstractMetric end

"""Prepared diagonal weights; construct through `FrequencyDomainMetric`."""
struct DiagonalMetric{W <: AbstractMatrix} <: AbstractMetric
    weights::W
end

"""
Prepared frequency blocks; `whiteners[:, :, k]` is the inverse lower Cholesky
factor of the one-sided spectral matrix at bin `k`. Construct through
`FrequencyDomainMetric`; no inverse spectral covariance matrix is formed.
"""
struct BlockDiagonalMetric{D, W <: AbstractArray{<:Number, 3}} <: AbstractMetric
    whiteners::W
end

BlockDiagonalMetric(w::AbstractArray{<:Number, 3}) = BlockDiagonalMetric{size(w, 1), typeof(w)}(w)

"""A prepared metric on an ordered set of channels and a uniform frequency grid."""
struct FrequencyDomainMetric{F, C, M <: AbstractMetric, T}
    frequencies::F
    channels::C
    metric::M
    active_bins::UnitRange{Int}
    normalization::T
end

function check_grid(f::AbstractRange)
    length(f) >= 2 || throw(ArgumentError("need at least two frequencies"))
    isfinite(first(f)) && isfinite(last(f)) && first(f) >= 0 && step(f) > 0 ||
        throw(ArgumentError("need a finite, nonnegative, increasing frequency range"))
    return nothing
end

function channel_order(channels)
    names = Tuple(channels)
    !isempty(names) && all(x -> x isa Symbol || x isa AbstractString, names) ||
        throw(ArgumentError("channels must be nonempty names (symbols or strings)"))
    canonical = map(Symbol, names)
    allunique(canonical) || throw(ArgumentError("channel names must be unique"))
    return canonical
end
