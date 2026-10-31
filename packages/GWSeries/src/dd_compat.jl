# DimensionalData compares numeric lookups with a relative `≈`, which accepts GPS epochs
# seconds apart. Replace the value comparison for our grids with the few-ulp grid check,
# which throws `ArgumentError` naming the quantity. Grids of different lengths are left to
# the generic method (singleton broadcasting, size errors). This hooks private DD API; the
# DimensionalData compat bound and the epoch-shift regression test guard it.
# DD 0.30 tuple comparison calls this hook directly, bypassing public pairwise dispatch.
for D in (:GPSTime, :Freq)
    @eval function Dimensions._comparedims2(a::$D, b::$D; val = false, msg = nothing, kw...)
        ga, gb = parent(lookup(a)), parent(lookup(b))
        if val && comparable_grids(ga, gb)
            if isnothing(msg)
                grids_agree(ga, gb) || return false
            else
                check_grid(ga, gb)
            end
            val = false
        end
        return invoke(Dimensions._comparedims2, Tuple{Dimension, Dimension}, a, b; val, msg, kw...)
    end
end

# DD 0.30 has no public terminal-index dispatch alias. Validate only after selector
# conversion so the comparison uses the selected destination block.
function Base.setindex!(s::AbstractSeries, v::AbstractDimArray, i1::DD.StandardIndices, I::DD.StandardIndices...)
    check_labelled_operands(view(s, i1, I...), v)
    setindex!(parent(s), parent(v), i1, I...)
    return s
end
# Resolve DD's overlapping DimIndices/empty-index-array selector method.
Base.setindex!(s::AbstractSeries, v::AbstractDimArray, i1::DD._DimIndicesAmb, I::DD._DimIndicesAmb...) =
    invoke(setindex!, Tuple{AbstractDimArray, Any, DD._DimIndicesAmb, Vararg{DD._DimIndicesAmb}}, s, v, i1, I...)
