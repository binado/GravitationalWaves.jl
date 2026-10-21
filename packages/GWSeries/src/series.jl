"""
    TimeSeries(array; epoch, channels, sampling_interval | sampling_frequency)

Rows are time samples and columns are ordered channels. Supply exactly one of
`sampling_interval` and `sampling_frequency`. The grid fields are immutable; only the
samples can change.
"""
struct TimeSeries{T, A <: AbstractMatrix{T}, E <: Real, D <: Real, C <: Tuple} <: AbstractSeries{T}
    data::A
    epoch::E
    sampling_interval::D
    channels::C
    function TimeSeries(a::AbstractMatrix{T}, m::NamedTuple, ::Validated) where {T}
        return new{T, typeof(a), typeof(m.epoch), typeof(m.sampling_interval), typeof(m.channels)}(
            a, m.epoch, m.sampling_interval, m.channels,
        )
    end
end

"""
    FrequencySeries(array; epoch, frequency_spacing, channels)

Complete one-sided spectrum with DC in the first row and ordered channels in columns.
Frequency coordinates are `(0:size(array, 1)-1) * frequency_spacing`; original time
length is not stored and must be supplied explicitly for inverse transforms.
"""
struct FrequencySeries{T, A <: AbstractMatrix{T}, E <: Real, D <: Real, C <: Tuple} <: AbstractSeries{T}
    data::A
    epoch::E
    frequency_spacing::D
    channels::C
    function FrequencySeries(a::AbstractMatrix{T}, m::NamedTuple, ::Validated) where {T}
        return new{T, typeof(a), typeof(m.epoch), typeof(m.frequency_spacing), typeof(m.channels)}(
            a, m.epoch, m.frequency_spacing, m.channels,
        )
    end
end

matrix_storage(a::AbstractMatrix{<:Number}) = a
matrix_storage(a::AbstractVector{<:Number}) = reshape(a, :, 1)
matrix_storage(a) = throw(ArgumentError("samples must be a numeric vector or matrix"))

function time_metadata(; epoch::Real, channels, sampling_interval = nothing, sampling_frequency = nothing)
    (sampling_interval === nothing) != (sampling_frequency === nothing) ||
        throw(ArgumentError("supply exactly one of sampling_interval and sampling_frequency"))
    t0 = metadata_value(epoch)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    spacing = if sampling_interval === nothing
        fs = metadata_value(sampling_frequency)
        isfinite(fs) && fs > 0 || throw(ArgumentError("sampling frequency must be finite and positive"))
        inv(fs)
    else
        metadata_value(sampling_interval)
    end
    isfinite(spacing) && spacing > 0 || throw(ArgumentError("sampling interval must be finite and positive"))
    return (; epoch = t0, sampling_interval = spacing, channels = channel_order(channels))
end

function TimeSeries(a; kwargs...)
    data = matrix_storage(a)
    m = time_metadata(; kwargs...)
    check_storage(data, m.channels)
    isfinite(m.epoch + max(size(data, 1) - 1, 0) * m.sampling_interval) ||
        throw(ArgumentError("sample coordinates must be finite"))
    return TimeSeries(data, m, Validated())
end

function frequency_metadata(; epoch::Real, channels, frequency_spacing::Real)
    t0, df = metadata_value(epoch), metadata_value(frequency_spacing)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    isfinite(df) && df > 0 || throw(ArgumentError("frequency spacing must be finite and positive"))
    return (; epoch = t0, frequency_spacing = df, channels = channel_order(channels))
end

function FrequencySeries(a; kwargs...)
    data = matrix_storage(a)
    m = frequency_metadata(; kwargs...)
    check_storage(data, m.channels)
    isfinite(max(size(data, 1) - 1, 0) * m.frequency_spacing) || throw(ArgumentError("frequency coordinates must be finite"))
    return FrequencySeries(data, m, Validated())
end

function check_storage(a, names)
    Base.require_one_based_indexing(a)
    size(a, 2) == length(names) || throw(DimensionMismatch("channel count does not match storage"))
    return nothing
end

"""Canonical constructor keywords, independent of sample storage and precision."""
metadata(s::TimeSeries) = (; epoch = s.epoch, sampling_interval = s.sampling_interval, channels = s.channels)
metadata(s::FrequencySeries) = (; epoch = s.epoch, frequency_spacing = s.frequency_spacing, channels = s.channels)
"""Return parent storage and constructor keywords without copying samples."""
flatten(s::AbstractSeries) = (parent(s), metadata(s))
Base.parent(s::AbstractSeries) = s.data
Base.size(s::AbstractSeries) = size(parent(s))
Base.axes(s::AbstractSeries) = axes(parent(s))
Base.IndexStyle(::Type{<:AbstractSeries}) = IndexCartesian()
Base.IndexStyle(::Type{<:TimeSeries{T, A}}) where {T, A} = IndexStyle(A)
Base.IndexStyle(::Type{<:FrequencySeries{T, A}}) where {T, A} = IndexStyle(A)
Base.@propagate_inbounds Base.getindex(s::AbstractSeries, i::Int) = parent(s)[i]
Base.@propagate_inbounds Base.getindex(s::AbstractSeries, i::Int, j::Int) = parent(s)[i, j]
Base.@propagate_inbounds function Base.setindex!(s::AbstractSeries, v, i::Int)
    parent(s)[i] = v
    return s
end
Base.@propagate_inbounds function Base.setindex!(s::AbstractSeries, v, i::Int, j::Int)
    parent(s)[i, j] = v
    return s
end
Base.getindex(s::AbstractSeries, i::CartesianIndex) = parent(s)[i]
Base.setindex!(s::AbstractSeries, v, i::CartesianIndex) = (setindex!(parent(s), v, i); s)
Base.dataids(s::AbstractSeries) = Base.dataids(parent(s))
Base.copy(s::AbstractSeries) = rebuild(s, copy(parent(s)), metadata(s))
# Base routines such as reverse, circshift and cumsum allocate through `similar` and then reorder or
# combine samples, so `similar` must not claim the source grid. Use `copy` for a series-valued buffer.
Base.similar(s::AbstractSeries, ::Type{T}) where {T} = similar(parent(s), T)
Base.similar(s::AbstractSeries, ::Type{T}, dims::Dims) where {T} = similar(parent(s), T, dims)
Base.similar(s::AbstractSeries, dims::Dims) = similar(parent(s), eltype(s), dims)
rebuild(::TimeSeries, a, m) = TimeSeries(a, m, Validated())
rebuild(::FrequencySeries, a, m) = FrequencySeries(a, m, Validated())

epoch(s::AbstractSeries) = s.epoch
channels(s::AbstractSeries) = s.channels
nsamples(s::AbstractSeries) = size(s, 1)
sampling_interval(s::TimeSeries) = s.sampling_interval
Base.step(s::TimeSeries) = sampling_interval(s)
sampling_frequency(s::TimeSeries) = inv(sampling_interval(s))
duration(s::TimeSeries) = nsamples(s) * sampling_interval(s)
times(s::TimeSeries) = range(epoch(s); step = sampling_interval(s), length = nsamples(s))
frequency_spacing(s::FrequencySeries) = s.frequency_spacing
frequencies(s::FrequencySeries) = range(zero(frequency_spacing(s)); step = frequency_spacing(s), length = nsamples(s))
coordinates(s::TimeSeries) = times(s)
coordinates(s::FrequencySeries) = frequencies(s)
