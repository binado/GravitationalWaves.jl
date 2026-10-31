module GWSeriesMakieExt

using GWSeries
import Makie

function representation_function(representation, ::Type{T}) where {T}
    representation === nothing && T <: Real && return identity
    representation === :magnitude && return abs
    representation === :phase && return angle
    representation === :real && return real
    representation === :imaginary && return imag
    throw(ArgumentError("specify representation as :magnitude, :phase, :real, or :imaginary for frequency data"))
end

plot_coordinates(s::TimeSeries) = range(zero(sampling_interval(s)); step = sampling_interval(s), length = size(s, 1))
plot_coordinates(s::FrequencySeries) = frequencies(s)
plot_transform(::TimeSeries, ::Nothing) = identity
plot_transform(s::FrequencySeries, representation) = representation_function(representation, eltype(s))
plot_xlabel(s::TimeSeries) = "Time since GPS $(epoch(s)) (s)"
plot_xlabel(::FrequencySeries) = "Frequency (Hz)"
plot_ylabel(::TimeSeries, _) = "Amplitude"
plot_ylabel(::FrequencySeries, ::Nothing) = "Amplitude"
plot_ylabel(::FrequencySeries, representation) = string(uppercasefirst(string(representation)))

function Makie.lines!(ax, s::Union{TimeSeries, FrequencySeries}; representation = nothing, kwargs...)
    transform = plot_transform(s, representation)
    isempty(channels(s)) && throw(ArgumentError("cannot plot a series without channels"))
    xs = plot_coordinates(s)
    first_plot = nothing
    for (i, channel) in pairs(channels(s))
        p = Makie.lines!(ax, xs, transform.(view(parent(s), :, i)); label = string(channel), kwargs...)
        isnothing(first_plot) && (first_plot = p)
    end
    return first_plot
end

function Makie.lines(s::Union{TimeSeries, FrequencySeries}; representation = nothing, figure = (;), axis = (;), kwargs...)
    # Validate before constructing a figure.
    plot_transform(s, representation)
    fig = Makie.Figure(; figure...)
    axis_options = merge((; xlabel = plot_xlabel(s), ylabel = plot_ylabel(s, representation)), axis)
    ax = Makie.Axis(fig[1, 1]; axis_options...)
    p = Makie.lines!(ax, s; representation, kwargs...)
    Makie.axislegend(ax)
    return Makie.FigureAxisPlot(fig, ax, p)
end

end
