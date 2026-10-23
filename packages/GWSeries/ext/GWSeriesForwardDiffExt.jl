module GWSeriesForwardDiffExt

using GWSeries, ForwardDiff

function GWSeries.physical_rfft(a::AbstractMatrix{ForwardDiff.Dual{Tag, V, N}}, dt::Real) where {Tag, V, N}
    values = GWSeries.physical_rfft(ForwardDiff.value.(a), dt)
    partials = ntuple(N) do k
        GWSeries.physical_rfft(map(x -> ForwardDiff.partials(x)[k], a), dt)
    end
    output = map(eachindex(values)) do i
        re = ForwardDiff.Dual{Tag}(real(values[i]), ntuple(k -> real(partials[k][i]), N))
        im = ForwardDiff.Dual{Tag}(imag(values[i]), ntuple(k -> imag(partials[k][i]), N))
        complex(re, im)
    end
    return reshape(output, size(values))
end

end
