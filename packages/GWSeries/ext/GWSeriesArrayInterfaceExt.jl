module GWSeriesArrayInterfaceExt
using GWSeries
import ArrayInterface
for S in (TimeSeries, FrequencySeries)
    @eval begin
        ArrayInterface.parent_type(::Type{<:$S{T, A}}) where {T, A} = A
        ArrayInterface.can_setindex(::Type{<:$S{T, A}}) where {T, A} = ArrayInterface.can_setindex(A)
        ArrayInterface.fast_scalar_indexing(::Type{<:$S{T, A}}) where {T, A} = ArrayInterface.fast_scalar_indexing(A)
        # Series do not promise a pointer interface, even when their parent has one.
        ArrayInterface.defines_strides(::Type{<:$S}) = false
    end
end
end
