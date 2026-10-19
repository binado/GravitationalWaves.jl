module GWTimeFrequencyDomainsFFTWExt

using GWTimeFrequencyDomains
using FFTW
using LinearAlgebra: mul!

const StridedIndices = Tuple{Vararg{Union{Int, AbstractRange{Int}}}}
const RealCPU = Union{Array{Float32}, Array{Float64}, SubArray{<:Union{Float32, Float64}, N, <:Array, <:StridedIndices} where {N}}
const ComplexCPU = Union{Array{ComplexF32}, Array{ComplexF64}, SubArray{<:Union{ComplexF32, ComplexF64}, N, <:Array, <:StridedIndices} where {N}}

function GWTimeFrequencyDomains.to_frequencydomain(d::TimeDomain, x::RealCPU)
    values = rfft(x, 1)
    values .*= sample_spacing(d)
    return values
end
function GWTimeFrequencyDomains.to_timedomain(d::FrequencyDomain, h::Union{RealCPU, ComplexCPU})
    n = 2 * (nsamples(d) - 1)
    return inverse_transform(h, n, n * frequency_spacing(d))
end
function GWTimeFrequencyDomains.to_timedomain(d::TimeDomain, h::Union{RealCPU, ComplexCPU})
    return inverse_transform(h, nsamples(d), sampling_frequency(d))
end
function inverse_transform(h, n, scale)
    values = irfft(h, n, 1)
    values .*= scale
    return values
end

# Inverse plans own a scratch input to guarantee source preservation, including
# FFTW's multidimensional/batched real-inverse implementation. Plans are not reentrant.
struct DomainFFTPlan{P, D, S, B, T, N}
    fft::P
    output::D
    scale::S
    buffer::B
    input_shape::NTuple{N, Int}
    output_shape::NTuple{N, Int}
end
function make_plan(fft, output, scale, buffer, ::Type{T}, input_shape::NTuple{N, Int}, output_shape::NTuple{N, Int}) where {T, N}
    return DomainFFTPlan{typeof(fft), typeof(output), typeof(scale), typeof(buffer), T, N}(fft, output, scale, buffer, input_shape, output_shape)
end
function GWTimeFrequencyDomains.plan_to_frequencydomain(d::TimeDomain, x::RealCPU)
    check_data(d, x)
    output = frequency_domain(d)
    shape = ntuple(i -> i == 1 ? nsamples(output) : size(x, i), ndims(x))
    return make_plan(plan_rfft(x, 1; flags = FFTW.ESTIMATE), output, sample_spacing(d), nothing, Complex{eltype(x)}, size(x), shape)
end
function GWTimeFrequencyDomains.plan_to_timedomain(d::FrequencyDomain, h::Union{RealCPU, ComplexCPU})
    check_data(d, h)
    output = time_domain(d)
    buffer = complex.(h)
    shape = ntuple(i -> i == 1 ? nsamples(output) : size(h, i), ndims(h))
    fft = plan_irfft(buffer, nsamples(output), 1; flags = FFTW.ESTIMATE)
    return make_plan(fft, output, nsamples(output) * frequency_spacing(d), buffer, real(eltype(h)), size(h), shape)
end
GWTimeFrequencyDomains.output_domain(p::DomainFFTPlan) = p.output
function GWTimeFrequencyDomains.transform(p::DomainFFTPlan{P, D, S, B, T}, input::Union{RealCPU, ComplexCPU}) where {P, D, S, B, T}
    destination = Array{T}(undef, p.output_shape)
    return transform!(destination, p, input)
end
function GWTimeFrequencyDomains.transform!(destination::Union{RealCPU, ComplexCPU}, p::DomainFFTPlan, input::Union{RealCPU, ComplexCPU})
    source = if p.buffer === nothing
        input
    else
        copyto!(p.buffer, input)
        p.buffer
    end
    mul!(destination, p.fft, source)
    destination .*= p.scale
    return destination
end

end
