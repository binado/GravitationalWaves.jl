module GWSeriesMooncakeExt

using GWSeries, AbstractFFTs, ChainRulesCore
using Mooncake: Mooncake, @is_primitive, MinimalCtx, CoDual, NoRData, primal, tangent

@is_primitive MinimalCtx Tuple{typeof(GWSeries.physical_rfft), Matrix{P}, P} where {P <: Union{Float32, Float64}}

function Mooncake.rrule!!(
        ::CoDual{typeof(GWSeries.physical_rfft)}, x::CoDual{Matrix{P}}, dt::CoDual{P},
    ) where {P <: Union{Float32, Float64}}
    a, scale = primal(x), primal(dt)
    raw, fft_pullback = ChainRulesCore.rrule(AbstractFFTs.rfft, a, 1)
    y = raw .* scale
    dy = zero(y)
    dx = tangent(x)
    function pullback(::NoRData)
        # Reuse AbstractFFTs' rule, including odd/even real-FFT endpoint weights.
        _, da, _ = fft_pullback(dy .* scale)
        dx .+= ChainRulesCore.unthunk(da)
        dscale = sum(real, conj.(dy) .* raw)
        return NoRData(), NoRData(), dscale
    end
    return CoDual(y, dy), pullback
end

end
