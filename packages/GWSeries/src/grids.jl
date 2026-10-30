"""
    TimeGrid(epoch, step, length)

Uniform absolute sample times, used as the [`GPSTime`](@ref) lookup of a
[`TimeSeries`](@ref). It wraps the Base range `range(epoch; step, length)`, returned by
`parent`, so its values are exactly those of that range. Slicing with a positive-step
integer range returns a `TimeGrid` wrapping the sliced range; any other index returns plain
values, so reordered or irregular selections stop being series. Coordinates are any
`Real`, so epochs that depend on differentiated parameters (dual numbers) are allowed.
"""
struct TimeGrid{T <: Real, R <: AbstractRange{T}} <: AbstractRange{T}
    times::R
    TimeGrid(times::AbstractRange{T}) where {T <: Real} = new{T, typeof(times)}(times)
end
function TimeGrid(epoch::Real, step::Real, len::Integer)
    t0, dt = promote(float(epoch), float(step))
    return TimeGrid(range(t0; step = dt, length = Int(len)))
end
Base.parent(g::TimeGrid) = g.times

"""
    SpectralGrid(spacing, length, epoch, time_samples)

Frequencies of a complete one-sided spectrum, used as the [`Freq`](@ref) lookup of a
[`FrequencySeries`](@ref). It wraps the Base range `range(0; step = spacing, length)`,
returned by `parent`, and also records the epoch that phases are referenced to and the
original time length (`nothing` when unknown). Only the full axis keeps these records:
partial bands return plain ranges.
"""
struct SpectralGrid{T <: Real, E <: Real, R <: AbstractRange{T}} <: AbstractRange{T}
    bins::R
    epoch::E
    # Zero marks an unrecorded length (recorded lengths are at least two). A `Nothing`-typed
    # field nested in the dims makes Mooncake's rule compilation segfault on Julia 1.12.
    recorded_samples::Int
    SpectralGrid(bins::AbstractRange{T}, epoch::E, recorded_samples::Int) where {T <: Real, E <: Real} =
        new{T, E, typeof(bins)}(bins, epoch, recorded_samples)
end
function SpectralGrid(spacing::Real, len::Integer, epoch::Real, time_samples::Union{Nothing, Integer})
    bins = range(zero(spacing); step = spacing, length = Int(len))
    return SpectralGrid(bins, epoch, Int(something(time_samples, 0)))
end
Base.parent(g::SpectralGrid) = g.bins
time_samples(g::SpectralGrid) = iszero(g.recorded_samples) ? nothing : g.recorded_samples

# Values come from the wrapped range, so grids agree bitwise with `parent(grid)`.
for G in (:TimeGrid, :SpectralGrid)
    @eval begin
        Base.first(g::$G) = first(parent(g))
        Base.step(g::$G) = step(parent(g))
        Base.last(g::$G) = last(parent(g))
        Base.length(g::$G) = length(parent(g))
        Base.size(g::$G) = size(parent(g))
        Base.getindex(g::$G, i::Integer) = parent(g)[i]
        Base.reverse(g::$G) = reverse(parent(g))
    end
end
function Base.getindex(g::TimeGrid, r::AbstractRange{<:Integer})
    times = parent(g)[r]
    step(r) > 0 && isfinite(first(times)) && isfinite(step(times)) || return times
    return TimeGrid(times)
end
function Base.getindex(g::SpectralGrid, r::AbstractRange{<:Integer})
    checkbounds(g, r)
    first(r) == 1 && step(r) == 1 && length(r) == length(g) && return g
    return parent(g)[r]
end

# Grids are immutable, so a view is the slice itself; this keeps `view(series, rows, :)` a series.
Base.view(g::Union{TimeGrid, SpectralGrid}, r::AbstractRange{<:Integer}) = g[r]

"""
    ChannelNames(names::NTuple{N,Symbol})

Ordered channel labels used as the [`Chan`](@ref) lookup. Every slice is again a
`ChannelNames`; uniqueness is checked when a series is rebuilt.
"""
struct ChannelNames{N} <: AbstractVector{Symbol}
    names::NTuple{N, Symbol}
end

Base.size(::ChannelNames{N}) where {N} = (N,)
Base.Tuple(c::ChannelNames) = c.names
Base.getindex(c::ChannelNames, i::Integer) = c.names[i]
Base.getindex(c::ChannelNames, ::Colon) = c
Base.getindex(c::ChannelNames, i::AbstractVector{<:Integer}) = ChannelNames(Tuple(collect(c.names)[i]))
Base.view(c::ChannelNames, i::Union{Colon, AbstractVector{<:Integer}}) = c[i]
