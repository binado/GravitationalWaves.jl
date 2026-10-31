using GWSeries
# Storage prototypes for tests that prepare grids independently of sample values.
time_prototype(; nsamples, channels, kwargs...) = TimeSeries(zeros(nsamples, length(channels)); channels, kwargs...)
function frequency_prototype(; frequencies, channels, kwargs...)
    first(frequencies) == 0 || throw(ArgumentError("a complete frequency prototype must start at DC"))
    return FrequencySeries(zeros(ComplexF64, length(frequencies), length(channels)); frequency_spacing = step(frequencies), channels, kwargs...)
end
