module GWSeriesFFTWExt

using GWSeries
using FFTW: rfft, irfft

# Views must be strided (integer and range indices); FFTW rejects index-vector views.
const StridedIndices = Tuple{Vararg{Union{Int, AbstractRange{Int}}}}
const RealCPU = Union{Array{Float32}, Array{Float64}, SubArray{<:Union{Float32, Float64}, N, <:Array, <:StridedIndices} where {N}}
const ComplexCPU = Union{Array{ComplexF32}, Array{ComplexF64}, SubArray{<:Union{ComplexF32, ComplexF64}, N, <:Array, <:StridedIndices} where {N}}

function GWSeries.to_frequencyseries(x::TimeSeries{<:RealCPU})
    n = nsamples(x)
    n >= 2 && iseven(n) || throw(ArgumentError("forward FFT needs even N >= 2"))
    values = rfft(parent(x), 1)
    values .*= sample_spacing(x)
    f = range(0.0; step = sampling_frequency(x) / n, length = n ÷ 2 + 1)
    return FrequencySeries(values; frequencies = f, epoch = epoch(x), channels = channels(x))
end

function GWSeries.to_timeseries(h::FrequencySeries{<:Union{RealCPU, ComplexCPU}})
    m = nsamples(h)
    m >= 2 && first(frequencies(h)) == 0 ||
        throw(ArgumentError("inverse FFT needs at least two bins starting at zero; restore omitted low bins explicitly"))
    n = 2 * (m - 1)
    fs = n * frequency_spacing(h)
    values = irfft(parent(h), n, 1)
    values .*= fs
    return TimeSeries(values; epoch = epoch(h), sampling_frequency = fs, channels = channels(h))
end

end
