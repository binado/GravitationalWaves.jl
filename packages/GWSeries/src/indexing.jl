# Resolve channel labels first; numerical indexing and views are delegated to storage.
channel_indices(s, j) = j
channel_indices(s, j::Tuple{}) = Int[]
channel_indices(s, j::Tuple{Vararg{Integer}}) = collect(j)
function channel_indices(s, j::Union{Symbol, AbstractString})
    i = findfirst(==(Symbol(j)), channels(s))
    i === nothing && throw(ArgumentError("unknown channel $j"))
    return i
end
channel_indices(s, j::Union{AbstractVector{<:Union{Symbol, AbstractString}}, Tuple{Vararg{Union{Symbol, AbstractString}}}}) = [channel_indices(s, label) for label in j]

function sliced(s, a, rows, cols)
    a isa AbstractMatrix || return a
    rows = rows isa Colon ? axes(s, 1) : rows
    rows isa AbstractRange{<:Integer} && step(rows) == 1 || return a
    cols = cols isa Colon ? axes(s, 2) : cols
    cols isa AbstractVector || return a
    indices = cols isa AbstractVector{Bool} ? findall(cols) : cols
    allunique(indices) && !isempty(indices) || return a
    names = Tuple(channels(s)[i] for i in indices)
    s isa FrequencySeries && rows != axes(s, 1) && return a
    origin = first(coordinates(s)) + (first(rows) - 1) * step(coordinates(s))
    spacing = step(rows) * step(coordinates(s))
    isfinite(origin) && isfinite(spacing) && spacing > 0 || throw(ArgumentError("sliced coordinates must have a finite origin and positive finite spacing"))
    m = if s isa TimeSeries
        merge(metadata(s), (; epoch = origin, sampling_interval = spacing, channels = names))
    else
        merge(metadata(s), (; channels = names))
    end
    return rebuild(s, a, m)
end

function Base.getindex(s::AbstractSeries, rows, cols)
    indices = channel_indices(s, cols)
    return sliced(s, getindex(parent(s), rows, indices), rows, indices)
end
function Base.view(s::AbstractSeries, rows, cols)
    indices = channel_indices(s, cols)
    return sliced(s, view(parent(s), rows, indices), rows, indices)
end
Base.getindex(s::AbstractSeries, i) = getindex(parent(s), i)
Base.view(s::AbstractSeries, i) = view(parent(s), i)
function Base.setindex!(s::AbstractSeries, v, rows, cols)
    if v isa AbstractSeries
        destination = s[rows, cols]
        destination isa AbstractSeries && check_same_grid(destination, v)
    end
    setindex!(parent(s), v isa AbstractSeries ? parent(v) : v, rows, channel_indices(s, cols))
    return s
end
Base.setindex!(s::AbstractSeries, v, i) = (setindex!(parent(s), v, i); s)
