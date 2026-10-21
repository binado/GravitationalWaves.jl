using GWSeries
# Storage prototypes for tests that prepare grids independently of sample values.
time_prototype(; nsamples, channels, kwargs...) = TimeSeries(zeros(nsamples, length(channels)); channels, kwargs...)
frequency_prototype(; frequencies, channels, kwargs...) = FrequencySeries(zeros(ComplexF64, length(frequencies), length(channels)); frequencies, channels, kwargs...)
