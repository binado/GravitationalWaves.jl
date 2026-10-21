struct SeriesStyle <: Base.Broadcast.AbstractArrayStyle{2} end
Base.BroadcastStyle(::Type{<:AbstractSeries}) = SeriesStyle()
Base.BroadcastStyle(::SeriesStyle, ::SeriesStyle) = SeriesStyle()
Base.BroadcastStyle(::SeriesStyle, ::Base.Broadcast.DefaultArrayStyle{N}) where {N} = SeriesStyle()
Base.BroadcastStyle(::Base.Broadcast.DefaultArrayStyle{N}, ::SeriesStyle) where {N} = SeriesStyle()
Base.BroadcastStyle(::SeriesStyle, ::Base.Broadcast.AbstractArrayStyle{N}) where {N} = SeriesStyle()
Base.BroadcastStyle(::Base.Broadcast.AbstractArrayStyle{N}, ::SeriesStyle) where {N} = SeriesStyle()

series_inputs(s::AbstractSeries) = (s,)
series_inputs(x) = ()
series_inputs(b::Base.Broadcast.Broadcasted) = series_inputs(b.args)
series_inputs(t::Tuple) = (series_inputs(first(t))..., series_inputs(Base.tail(t))...)
series_inputs(::Tuple{}) = ()
function validate_inputs(inputs::Tuple)
    isempty(inputs) && return nothing
    for i in eachindex(inputs)
        any(j -> inputs[j] === inputs[i], 1:(i - 1)) && continue
        check_same_grid(first(inputs), inputs[i])
    end
    return nothing
end
unwrap(s::AbstractSeries) = parent(s)
unwrap(x) = x
unwrap(b::Base.Broadcast.Broadcasted) = Base.Broadcast.broadcasted(b.f, map(unwrap, b.args)...)
wrap_result(s::AbstractSeries, a) = a isa AbstractMatrix && size(a) == size(s) ? rebuild(s, a, metadata(s)) : a
function Base.copy(b::Base.Broadcast.Broadcasted{SeriesStyle})
    inputs = series_inputs(b)
    validate_inputs(inputs)
    return wrap_result(first(inputs), Base.copy(Base.Broadcast.instantiate(unwrap(b))))
end
function series_copyto!(dest, b)
    inputs = dest isa AbstractSeries ? (dest, series_inputs(b)...) : series_inputs(b)
    validate_inputs(inputs)
    copyto!(unwrap(dest), Base.Broadcast.instantiate(unwrap(b)))
    return dest
end
Base.copyto!(dest::AbstractArray, b::Base.Broadcast.Broadcasted{SeriesStyle}) = series_copyto!(dest, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted{SeriesStyle}) = series_copyto!(dest, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted) = series_copyto!(dest, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted{Nothing}) = series_copyto!(dest, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted{<:Base.Broadcast.AbstractArrayStyle{0}}) = series_copyto!(dest, b)

function array_operation(f, a, b)
    inputs = (series_inputs(a)..., series_inputs(b)...)
    validate_inputs(inputs)
    return wrap_result(first(inputs), f(unwrap(a), unwrap(b)))
end
for op in (:+, :-, :*)
    @eval begin
        Base.$op(a::AbstractSeries, b::AbstractMatrix) = array_operation($op, a, b)
        Base.$op(a::AbstractMatrix, b::AbstractSeries) = array_operation($op, a, b)
        Base.$op(a::AbstractSeries, b::AbstractSeries) = array_operation($op, a, b)
    end
end
Base.:*(a::AbstractSeries, b::AbstractVector) = parent(a) * b
Base.:*(a::AbstractSeries, b::Number) = wrap_result(a, parent(a) * b)
Base.:*(a::Number, b::AbstractSeries) = wrap_result(b, a * parent(b))
Base.:/(a::AbstractSeries, b::Number) = wrap_result(a, parent(a) / b)
Base.:-(a::AbstractSeries) = wrap_result(a, -parent(a))
Base.:+(a::AbstractSeries) = copy(a)
Base.map(f, a::AbstractSeries) = wrap_result(a, map(f, parent(a)))
Base.map(f, a::AbstractSeries, b::AbstractArray) = array_operation((x, y) -> map(f, x, y), a, b)
Base.map(f, a::AbstractArray, b::AbstractSeries) = array_operation((x, y) -> map(f, x, y), a, b)
Base.map(f, a::AbstractSeries, b::AbstractSeries) = array_operation((x, y) -> map(f, x, y), a, b)
Base.reshape(s::AbstractSeries, dims::Union{Int, Colon}...) = reshape(parent(s), dims...)
Base.reshape(s::AbstractSeries, dims::Dims) = reshape(parent(s), dims)
Base.reshape(s::AbstractSeries, dims::Tuple{Integer, Vararg{Integer}}) = reshape(parent(s), dims)
Base.reshape(s::AbstractSeries, dims::Tuple{Vararg{Union{Colon, Integer}}}) = reshape(parent(s), dims)
Base.reshape(s::AbstractSeries, dims::Tuple{Union{Integer, Base.OneTo}, Vararg{Union{Integer, Base.OneTo}}}) = reshape(parent(s), dims)
Base.vec(s::AbstractSeries) = vec(parent(s))
Base.transpose(s::AbstractSeries) = transpose(parent(s))
Base.adjoint(s::AbstractSeries) = adjoint(parent(s))
Base.permutedims(s::AbstractSeries, perm = (2, 1)) = permutedims(parent(s), perm)
Base.sum(s::AbstractSeries; kwargs...) = sum(parent(s); kwargs...)
Base.sum(f, s::AbstractSeries; kwargs...) = sum(f, parent(s); kwargs...)
Base.prod(s::AbstractSeries; kwargs...) = prod(parent(s); kwargs...)
Base.maximum(s::AbstractSeries; kwargs...) = maximum(parent(s); kwargs...)
Base.minimum(s::AbstractSeries; kwargs...) = minimum(parent(s); kwargs...)

function series_map!(f, destination, inputs...)
    series = (series_inputs(destination)..., series_inputs(inputs)...)
    validate_inputs(series)
    map!(f, unwrap(destination), map(unwrap, inputs)...)
    return destination
end
Base.map!(f, destination::AbstractSeries, inputs::AbstractArray...) = series_map!(f, destination, inputs...)
Base.map!(f, destination::AbstractArray, input::AbstractSeries, inputs::AbstractArray...) = series_map!(f, destination, input, inputs...)
Base.map!(f, destination::AbstractSeries, input::AbstractSeries, inputs::AbstractArray...) = series_map!(f, destination, input, inputs...)
function Base.copyto!(destination::AbstractSeries, source::AbstractSeries)
    check_same_grid(destination, source)
    copyto!(parent(destination), parent(source))
    return destination
end

Base.:*(a::Transpose{<:Any, <:AbstractVector}, b::AbstractSeries) = array_operation(*, a, b)
Base.:*(a::Adjoint{<:Any, <:AbstractVector}, b::AbstractSeries) = array_operation(*, a, b)
Base.prod(f, s::AbstractSeries; kwargs...) = prod(f, parent(s); kwargs...)
Base.maximum(f, s::AbstractSeries; kwargs...) = maximum(f, parent(s); kwargs...)
Base.minimum(f, s::AbstractSeries; kwargs...) = minimum(f, parent(s); kwargs...)
Base.extrema(s::AbstractSeries; kwargs...) = extrema(parent(s); kwargs...)
Base.extrema(f, s::AbstractSeries; kwargs...) = extrema(f, parent(s); kwargs...)
