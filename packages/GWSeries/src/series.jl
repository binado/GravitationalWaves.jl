const TimeDims = Tuple{GPSTime{<:Sampled{<:Any, <:TimeGrid}}, Chan{<:Categorical{Symbol, <:ChannelNames}}}
const FrequencyDims = Tuple{Freq{<:Sampled{<:Any, <:SpectralGrid}}, Chan{<:Categorical{Symbol, <:ChannelNames}}}

function check_dims(a::AbstractMatrix, grid, names::ChannelNames)
    Base.require_one_based_indexing(a)
    size(a) == (length(grid), length(names)) || throw(DimensionMismatch("storage size $(size(a)) does not match the grid and channels"))
    allunique(names) || throw(ArgumentError("channel names must be unique"))
    return nothing
end

"""
    TimeSeries(array; epoch, channels, sampling_interval | sampling_frequency)

Rows are time samples on the absolute [`GPSTime`](@ref) grid and columns are ordered
[`Chan`](@ref) channels. Supply exactly one of `sampling_interval` and
`sampling_frequency`. The grid is immutable; only the samples can change.
"""
struct TimeSeries{T, A <: AbstractMatrix{T}, D <: TimeDims} <: AbstractSeries{T, D, A}
    data::A
    dims::D
    function TimeSeries(a::AbstractMatrix{T}, grid::TimeGrid, names::ChannelNames) where {T}
        check_dims(a, grid, names)
        d = (GPSTime(grid_lookup(grid)), Chan(channel_lookup(names)))
        return new{T, typeof(a), typeof(d)}(a, d)
    end
end

"""
    FrequencySeries(array; epoch, frequency_spacing, channels)

Complete one-sided spectrum with DC in the first row and ordered channels in columns.
Frequency coordinates are `(0:size(array, 1)-1) * frequency_spacing`; original time
length is recorded in `time_samples` when the spectrum comes from a real FFT (it fixes
the even/odd ambiguity of the inverse) and is `nothing` otherwise, in which case inverse
transforms need an explicit length.
"""
struct FrequencySeries{T, A <: AbstractMatrix{T}, D <: FrequencyDims} <: AbstractSeries{T, D, A}
    data::A
    dims::D
    function FrequencySeries(a::AbstractMatrix{T}, grid::SpectralGrid, names::ChannelNames) where {T}
        check_dims(a, grid, names)
        d = (Freq(grid_lookup(grid)), Chan(channel_lookup(names)))
        return new{T, typeof(a), typeof(d)}(a, d)
    end
end

TimeSeries(a::AbstractMatrix, m::NamedTuple, ::Validated) =
    TimeSeries(a, TimeGrid(m.epoch, m.sampling_interval, size(a, 1)), ChannelNames(m.channels))
FrequencySeries(a::AbstractMatrix, m::NamedTuple, ::Validated) =
    FrequencySeries(a, SpectralGrid(m.frequency_spacing, size(a, 1), m.epoch, m.time_samples), ChannelNames(m.channels))
"""Rebuild a series of the same kind from storage and validated [`grid_metadata`](@ref)."""
reconstruct(::TimeSeries, a, m) = TimeSeries(a, m, Validated())
reconstruct(::FrequencySeries, a, m) = FrequencySeries(a, m, Validated())

matrix_storage(a::AbstractMatrix{<:Number}) = a
matrix_storage(a::AbstractVector{<:Number}) = reshape(a, :, 1)
matrix_storage(a) = throw(ArgumentError("samples must be a numeric vector or matrix"))

"""Validated `(epoch, sampling_interval, channels)` keywords; exactly one of `sampling_interval` and `sampling_frequency` is required."""
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

"""Validated `(epoch, frequency_spacing, channels, time_samples)` keywords for a [`FrequencySeries`](@ref)."""
function frequency_metadata(; epoch::Real, channels, frequency_spacing::Real, time_samples::Union{Nothing, Integer} = nothing)
    t0, df = metadata_value(epoch), metadata_value(frequency_spacing)
    isfinite(t0) || throw(ArgumentError("epoch must be finite"))
    isfinite(df) && df > 0 || throw(ArgumentError("frequency spacing must be finite and positive"))
    time_samples === nothing || time_samples >= 2 || throw(ArgumentError("time_samples must be at least two"))
    return (; epoch = t0, frequency_spacing = df, channels = channel_order(channels), time_samples = time_samples === nothing ? nothing : Int(time_samples))
end

function FrequencySeries(a; kwargs...)
    data = matrix_storage(a)
    m = frequency_metadata(; kwargs...)
    check_storage(data, m.channels)
    m.time_samples === nothing || size(data, 1) == m.time_samples ÷ 2 + 1 ||
        throw(DimensionMismatch("a real FFT of $(m.time_samples) samples has $(m.time_samples ÷ 2 + 1) bins, got $(size(data, 1))"))
    isfinite(max(size(data, 1) - 1, 0) * m.frequency_spacing) || throw(ArgumentError("frequency coordinates must be finite"))
    return FrequencySeries(data, m, Validated())
end

function check_storage(a, names)
    Base.require_one_based_indexing(a)
    size(a, 2) == length(names) || throw(DimensionMismatch("channel count does not match storage"))
    return nothing
end

# DimensionalData interface. Series carry no name, refdims, or DD metadata.
DD.refdims(::AbstractSeries) = ()
DD.name(::AbstractSeries) = NoName()
DD.metadata(::AbstractSeries) = NoMetadata()

# Every DD result path (broadcast, similar, slicing, reductions, products) rebuilds from the
# source array, so this one method decides whether a result is still a series.
DD.rebuild(::AbstractSeries, data, dims::Tuple, refdims, name, metadata) =
    rebuild_series(data, dims, refdims, name, metadata)
DD.rebuild(A::AbstractSeries; data = parent(A), dims = DD.dims(A), refdims = (), name = NoName(), metadata = NoMetadata(), kw...) =
    rebuild_series(data, dims, refdims, name, metadata)
rebuild_series(data::AbstractMatrix, dims::TimeDims, refdims, name, metadata) =
    allunique(parent(lookup(dims[2]))) ? TimeSeries(data, parent(lookup(dims[1])), parent(lookup(dims[2]))) :
    DimArray(data, dims, refdims, name, metadata)
rebuild_series(data::AbstractMatrix, dims::FrequencyDims, refdims, name, metadata) =
    allunique(parent(lookup(dims[2]))) ? FrequencySeries(data, parent(lookup(dims[1])), parent(lookup(dims[2]))) :
    DimArray(data, dims, refdims, name, metadata)
rebuild_series(data, dims, refdims, name, metadata) = DimArray(data, dims, refdims, name, metadata)

"""Drop the series guarantees, keeping samples and labels."""
dimarray(s::AbstractSeries) = DimArray(parent(s), dims(s))
# DimensionalData keeps the source dims for these, which would mislabel reordered samples.
for S in (:Real, :(Tuple{Vararg{Integer}}), :Any)
    @eval Base.circshift(s::AbstractSeries, shifts::$S) = DimArray(circshift(parent(s), shifts), dims(s))
end
Base.diff(s::AbstractSeries; dims) = dimarray(invoke(diff, Tuple{AbstractDimArray}, s; dims))
# Accumulating across channels mixes them; only accumulation along samples stays a series.
function Base.cumsum(s::AbstractSeries; dims)
    out = invoke(cumsum, Tuple{AbstractDimArray}, s; dims)
    return dimnum(s, dims) == 1 ? out : dimarray(out)
end
# DimensionalData's `map` skips lookup values; compare grids as broadcasting does.
function Base.map(f, A1::AbstractSeries, As::AbstractDimArray...)
    DD.comparedims(A1, As...; val = true)
    return invoke(map, Tuple{Any, AbstractDimArray, Vararg{AbstractDimArray}}, f, A1, As...)
end
# Assigning a series checks the destination block's grid; selectors reach this after conversion.
function Base.setindex!(s::AbstractSeries, v::AbstractSeries, i1::DD.StandardIndices, I::DD.StandardIndices...)
    destination = view(s, i1, I...)
    if destination isa AbstractSeries
        check_same_grid(destination, v)
    else
        DD.comparedims(destination, v; val = true)
    end
    setindex!(parent(s), parent(v), i1, I...)
    return s
end
# Resolves an ambiguity with DimensionalData's method for DimIndices and empty index arrays.
Base.setindex!(s::AbstractSeries, v::AbstractSeries, i1::DD._DimIndicesAmb, I::DD._DimIndicesAmb...) =
    invoke(setindex!, Tuple{AbstractDimArray, Any, DD._DimIndicesAmb, Vararg{DD._DimIndicesAmb}}, s, v, i1, I...)
function Base.copyto!(destination::AbstractSeries, source::AbstractSeries)
    check_same_grid(destination, source)
    copyto!(parent(destination), parent(source))
    return destination
end

"""Canonical constructor keywords, independent of sample storage and precision."""
grid_metadata(s::TimeSeries) = (; epoch = epoch(s), sampling_interval = sampling_interval(s), channels = channels(s))
grid_metadata(s::FrequencySeries) = (; epoch = epoch(s), frequency_spacing = frequency_spacing(s), channels = channels(s), time_samples = time_samples(s))
"""Return parent storage and constructor keywords without copying samples."""
flatten(s::AbstractSeries) = (parent(s), grid_metadata(s))

# The grids stay inside the lookups; public accessors return their raw Base ranges, so code
# consuming coordinates gets exactly the range arithmetic it had before DimensionalData.
grid(s::AbstractSeries) = parent(lookup(s, 1))
"""Absolute sample times as a Base range (`parent` of the [`TimeGrid`](@ref) lookup)."""
times(s::TimeSeries) = parent(grid(s))
"""Bin frequencies as a Base range (`parent` of the [`SpectralGrid`](@ref) lookup)."""
frequencies(s::FrequencySeries) = parent(grid(s))
coordinates(s::TimeSeries) = times(s)
coordinates(s::FrequencySeries) = frequencies(s)
epoch(s::TimeSeries) = first(times(s))
epoch(s::FrequencySeries) = grid(s).epoch
channels(s::AbstractSeries) = Tuple(parent(lookup(s, 2)))
nsamples(s::AbstractSeries) = size(s, 1)
sampling_interval(s::TimeSeries) = step(times(s))
Base.step(s::TimeSeries) = sampling_interval(s)
sampling_frequency(s::TimeSeries) = inv(sampling_interval(s))
duration(s::TimeSeries) = nsamples(s) * sampling_interval(s)
frequency_spacing(s::FrequencySeries) = step(frequencies(s))
"""Original time length of a spectrum produced by a real FFT, or `nothing` when unrecorded."""
time_samples(s::FrequencySeries) = time_samples(grid(s))
