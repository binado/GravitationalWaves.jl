"""
    TimeGrid(epoch, step, length)

Uniform absolute sample times `epoch + (i - 1) * step`, used as the [`GPSTime`](@ref)
lookup of a [`TimeSeries`](@ref). Slicing with a positive-step integer range returns a
`TimeGrid` with an advanced epoch and a multiplied step; any other index returns plain
values, so reordered or irregular selections stop being series. Coordinates are any
`Real`, so epochs that depend on differentiated parameters (dual numbers) are allowed.
"""
struct TimeGrid{T <: Real} <: AbstractRange{T}
    epoch::T
    step::T
    len::Int
end
TimeGrid(epoch::Real, step::Real, len::Integer) = TimeGrid(promote(float(epoch), float(step))..., Int(len))

Base.first(g::TimeGrid) = g.epoch
Base.step(g::TimeGrid) = g.step
Base.length(g::TimeGrid) = g.len
Base.size(g::TimeGrid) = (g.len,)
Base.last(g::TimeGrid) = sample_time(g, g.len)
# One rounding of `epoch + k * step`, as in the double-double `StepRangeLen` it replaces.
sample_time(g::TimeGrid, i::Integer) = fma(oftype(g.step, i - 1), g.step, g.epoch)
function Base.getindex(g::TimeGrid, i::Integer)
    checkbounds(g, i)
    return sample_time(g, i)
end
function Base.getindex(g::TimeGrid, r::AbstractRange{<:Integer})
    checkbounds(g, r)
    step(r) > 0 || return [sample_time(g, i) for i in r]
    origin, spacing = sample_time(g, first(r)), step(r) * step(g)
    isfinite(origin) && isfinite(spacing) || return [sample_time(g, i) for i in r]
    return TimeGrid(origin, spacing, length(r))
end
Base.reverse(g::TimeGrid) = StepRangeLen(last(g), -step(g), length(g))

"""
    SpectralGrid(spacing, length, epoch, time_samples)

Frequencies `(i - 1) * spacing` of a complete one-sided spectrum, used as the
[`Freq`](@ref) lookup of a [`FrequencySeries`](@ref). It also records the epoch that
phases are referenced to and the original time length (`nothing` when unknown). Only the
full axis keeps these records: partial bands return plain ranges.
"""
struct SpectralGrid{T <: Real, E <: Real} <: AbstractRange{T}
    spacing::T
    len::Int
    epoch::E
    # Zero marks an unrecorded length (recorded lengths are at least two). A `Nothing`-typed
    # field nested in the dims makes Mooncake's rule compilation segfault on Julia 1.12.
    recorded_samples::Int
end
SpectralGrid(spacing::Real, len::Integer, epoch::Real, time_samples::Union{Nothing, Integer}) =
    SpectralGrid(spacing, Int(len), epoch, something(time_samples, 0))
time_samples(g::SpectralGrid) = iszero(g.recorded_samples) ? nothing : g.recorded_samples

Base.first(g::SpectralGrid) = zero(g.spacing)
Base.step(g::SpectralGrid) = g.spacing
Base.length(g::SpectralGrid) = g.len
Base.size(g::SpectralGrid) = (g.len,)
Base.last(g::SpectralGrid) = bin_frequency(g, g.len)
bin_frequency(g::SpectralGrid, i::Integer) = oftype(g.spacing, i - 1) * g.spacing
function Base.getindex(g::SpectralGrid, i::Integer)
    checkbounds(g, i)
    return bin_frequency(g, i)
end
function Base.getindex(g::SpectralGrid, r::AbstractRange{<:Integer})
    checkbounds(g, r)
    first(r) == 1 && step(r) == 1 && length(r) == length(g) && return g
    return StepRangeLen(bin_frequency(g, first(r)), step(r) * g.spacing, length(r))
end
Base.reverse(g::SpectralGrid) = StepRangeLen(last(g), -step(g), length(g))

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
