"""Absolute GPS sample-time dimension of a [`TimeSeries`](@ref); select with `x[GPSTime = t0 .. t1]`."""
@dim GPSTime TimeDim "GPS time"
"""Frequency dimension of a [`FrequencySeries`](@ref); select with `h[Freq = 20 .. 1024]`."""
@dim Freq "Frequency"
"""Channel dimension of a series; select with `x[Chan = At(:H1)]` or `x[Chan = At([:L1, :H1])]`."""
@dim Chan "Channel"

grid_lookup(g::Union{TimeGrid, SpectralGrid}) =
    Sampled(g; order = ForwardOrdered(), span = Regular(step(g)), sampling = Points(), metadata = NoMetadata())
# Unordered: channel order is a label, not a sort order, so swapped channels report a value mismatch.
channel_lookup(c::ChannelNames) = Categorical(c; order = Unordered(), metadata = NoMetadata())

check_grid(a::TimeGrid, b::TimeGrid) = check_range(a, b, "epochs", EPOCH_HINT)
function check_grid(a::SpectralGrid, b::SpectralGrid)
    check_range(a, b, "origins")
    ta, tb = time_samples(a), time_samples(b)
    ta === nothing || tb === nothing || ta == tb || throw(ArgumentError("original time lengths differ: $ta vs $tb"))
    epoch_spacing = inv(something(ta, tb, max(2 * (length(a) - 1), 1)) * step(a))
    check_coordinate(a.epoch, b.epoch, epoch_spacing, "epochs", EPOCH_HINT)
    return nothing
end

# Only equal-length grids of one kind are compared here; everything else takes the generic path.
comparable_grids(a::TimeGrid, b::TimeGrid) = length(a) == length(b)
comparable_grids(a::SpectralGrid, b::SpectralGrid) = length(a) == length(b)
comparable_grids(a, b) = false

function grids_agree(a, b)
    try
        check_grid(a, b)
    catch e
        e isa Union{ArgumentError, DimensionMismatch} || rethrow()
        return false
    end
    return true
end

# DimensionalData compares numeric lookups with a relative `≈`, which accepts GPS epochs
# seconds apart. Replace the value comparison for our grids with the few-ulp grid check,
# which throws `ArgumentError` naming the quantity. Grids of different lengths are left to
# the generic method (singleton broadcasting, size errors). This hooks private DD API; the
# DimensionalData compat bound and the epoch-shift regression test guard it.
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
