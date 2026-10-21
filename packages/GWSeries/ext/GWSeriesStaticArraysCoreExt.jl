module GWSeriesStaticArraysCoreExt
using GWSeries
using StaticArraysCore: StaticArray, StaticArrayStyle

for op in (:+, :-)
    @eval begin
        Base.$op(a::StaticArray{<:Tuple, <:Any, 2}, b::AbstractSeries) = GWSeries.array_operation($op, a, b)
        Base.$op(a::AbstractSeries, b::StaticArray{<:Tuple, <:Any, 2}) = GWSeries.array_operation($op, a, b)
    end
end
Base.map(f, a::StaticArray, b::AbstractSeries) = GWSeries.array_operation((x, y) -> map(f, x, y), a, b)
Base.map(f, a::AbstractSeries, b::StaticArray) = GWSeries.array_operation((x, y) -> map(f, x, y), a, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted{<:StaticArrayStyle}) = GWSeries.series_copyto!(dest, b)
Base.copyto!(dest::AbstractSeries, b::Base.Broadcast.Broadcasted{StaticArrayStyle{0}}) = GWSeries.series_copyto!(dest, b)
end
