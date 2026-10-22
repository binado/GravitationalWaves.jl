function check_fourier_domain(s::TimeSeries)
    nsamples(s) >= 2 || throw(ArgumentError("FFT needs at least two time samples"))
    return nothing
end
function check_fourier_domain(s::FrequencySeries; nsamples::Integer)
    nsamples >= 2 && size(s, 1) == nsamples ÷ 2 + 1 ||
        throw(ArgumentError("inverse real FFT needs N÷2+1 bins for its output length"))
    return nothing
end

# AbstractFFTs only defines the interface; without a backend its generic fallbacks recurse
# (a StackOverflowError for rfft) or raise an unhelpful MethodError. Detect that on the host.
function backend_loaded(f, types::Type{<:Tuple})
    hasmethod(f, types) || return false
    return parentmodule(which(f, types)) !== AbstractFFTs
end
function require_backend(a::AbstractArray, inverse::Bool)
    loaded = if inverse
        backend_loaded(AbstractFFTs.plan_brfft, Tuple{typeof(a), Int, Int})
    else
        backend_loaded(AbstractFFTs.plan_rfft, Tuple{typeof(a), Int})
    end
    loaded || throw(ArgumentError("no AbstractFFTs backend is loaded for $(typeof(a)); load one, for example `using FFTW` for CPU arrays"))
    return nothing
end

function transform_dims(dims)
    (dims == 1 || dims == (1,) || dims == [1]) || throw(ArgumentError("series FFTs transform only the sample dimension (1)"))
    return 1
end

function fourier_metadata(s::TimeSeries)
    return (; epoch = epoch(s), frequency_spacing = inv(duration(s)), channels = channels(s), time_samples = nsamples(s))
end

# The recorded length fixes the even/odd ambiguity of the inverse (N and N+1 can share a bin count).
function inverse_length(s::FrequencySeries, n)
    known = time_samples(s)
    n === nothing && known === nothing && throw(ArgumentError("the original time length is not recorded; pass it explicitly"))
    n === nothing && return known
    known === nothing || known == n || throw(ArgumentError("requested time length $n differs from the recorded length $known"))
    return Int(n)
end
function inverse_metadata(s::FrequencySeries, n)
    return (; epoch = epoch(s), sampling_interval = inv(n * frequency_spacing(s)), channels = channels(s))
end

# Array boundary for optional AD extensions. Metadata stays outside the primitive.
function physical_rfft(a::AbstractMatrix{<:Real}, dt::Real)
    require_backend(a, false)
    return AbstractFFTs.rfft(a, 1) .* dt
end

"""Physical dt*rfft along rows; the loaded AbstractFFTs backend owns execution."""
function AbstractFFTs.rfft(s::TimeSeries{<:Real}, dims = 1)
    transform_dims(dims)
    check_fourier_domain(s)
    a = physical_rfft(parent(s), sampling_interval(s))
    return FrequencySeries(a, fourier_metadata(s), Validated())
end
AbstractFFTs.rfft(s::TimeSeries, dims = 1) = throw(ArgumentError("wrapped forward transforms require real samples"))

AbstractFFTs.irfft(s::FrequencySeries) = AbstractFFTs.irfft(s, inverse_length(s, nothing))
function AbstractFFTs.irfft(s::FrequencySeries, n::Integer, dims = 1)
    transform_dims(dims)
    n = inverse_length(s, n)
    check_fourier_domain(s; nsamples = n)
    spectrum = complex.(parent(s))
    require_backend(spectrum, true)
    a = AbstractFFTs.irfft(spectrum, n, 1)
    a .*= n * frequency_spacing(s)
    return TimeSeries(a, inverse_metadata(s, n), Validated())
end
for f in (:fft, :ifft, :bfft, :plan_fft, :plan_ifft, :plan_bfft)
    @eval AbstractFFTs.$f(s::AbstractSeries, args...; kwargs...) = throw(ArgumentError("use parent(series) for complex or two-sided transforms"))
end
for f in (:fft!, :ifft!, :bfft!, :plan_fft!, :plan_ifft!, :plan_bfft!)
    @eval AbstractFFTs.$f(s::AbstractSeries, args...; kwargs...) = throw(ArgumentError("in-place FFTs cannot change a series type; use a physical transform plan or parent(series)"))
end
for f in (:fftshift, :ifftshift)
    @eval AbstractFFTs.$f(s::AbstractSeries, args...) = throw(ArgumentError("shifting bins changes series coordinates; use parent(series) or the physical FFT methods"))
end

"""Transform a time series with physical normalization (real FFT for real samples)."""
to_frequencydomain(s::TimeSeries{<:Real}) = AbstractFFTs.rfft(s)
to_frequencydomain(s::TimeSeries) = throw(ArgumentError("wrapped forward transforms require real samples"))
"""Invert a physical one-sided spectrum; `nsamples` is required unless the spectrum records its time length."""
to_timedomain(s::FrequencySeries; nsamples::Union{Nothing, Integer} = nothing) = AbstractFFTs.irfft(s, inverse_length(s, nsamples))

struct Forward end
struct Inverse end

"""Reusable AbstractFFTs plan and private scratch storage; not reentrant."""
struct SeriesFFTPlan{K <: Union{Forward, Inverse}, P, I, W, IM, OM, S}
    fft::P
    input::I
    work::W
    input_metadata::IM
    output_metadata::OM
    scale::S
end
function SeriesFFTPlan{K}(fft::P, input::I, work::W, input_metadata::IM, output_metadata::OM, scale::S) where {K, P, I, W, IM, OM, S}
    return SeriesFFTPlan{K, P, I, W, IM, OM, S}(fft, input, work, input_metadata, output_metadata, scale)
end
metadata(p::SeriesFFTPlan) = p.output_metadata
input_type(::SeriesFFTPlan{Forward}) = TimeSeries
input_type(::SeriesFFTPlan{Inverse}) = FrequencySeries
output_type(::SeriesFFTPlan{Forward}) = FrequencySeries
output_type(::SeriesFFTPlan{Inverse}) = TimeSeries

"""Prepare a physical forward FFT. Backend-specific keywords pass to AbstractFFTs."""
function plan_to_frequencydomain(s::TimeSeries{<:Real}; kwargs...)
    check_fourier_domain(s)
    input = copy(parent(s))
    require_backend(input, false)
    fft = AbstractFFTs.plan_rfft(input, 1; kwargs...)
    shape = (nsamples(s) ÷ 2 + 1, size(s, 2))
    work = similar(input, Complex{real(eltype(input))}, shape)
    return SeriesFFTPlan{Forward}(fft, input, work, metadata(s), fourier_metadata(s), sampling_interval(s))
end
plan_to_frequencydomain(s::TimeSeries; kwargs...) = throw(ArgumentError("wrapped forward transforms require real samples"))

"""Prepare a source-preserving inverse real FFT; `nsamples` is required unless the spectrum records its time length."""
function plan_to_timedomain(s::FrequencySeries; nsamples::Union{Nothing, Integer} = nothing, kwargs...)
    nsamples = inverse_length(s, nsamples)
    check_fourier_domain(s; nsamples)
    input = complex.(parent(s))
    require_backend(input, true)
    fft = AbstractFFTs.plan_irfft(input, nsamples, 1; kwargs...)
    work = similar(input, real(eltype(input)), (nsamples, size(s, 2)))
    return SeriesFFTPlan{Inverse}(fft, input, work, metadata(s), inverse_metadata(s, nsamples), nsamples * frequency_spacing(s))
end

function transform_storage!(destination::AbstractMatrix, p::SeriesFFTPlan, input::AbstractMatrix)
    Base.mightalias(destination, input) && throw(ArgumentError("transform destination must not alias its input"))
    size(input) == size(p.input) && size(destination) == size(p.work) || throw(DimensionMismatch("transform storage does not match the prepared shapes"))
    copyto!(p.input, input)
    mul!(p.work, p.fft, p.input)
    destination .= p.work .* p.scale
    return destination
end
transform(p::SeriesFFTPlan, input::AbstractMatrix) = transform_storage!(similar(p.work), p, input)

# The plan's metadata is already validated, so compare grids through the trusted constructor.
function check_plan_series(::Type{S}, s, meta, role) where {S}
    s isa S || throw(ArgumentError("series kind does not match transform $role"))
    check_same_grid(rebuild(s, parent(s), meta), s)
    return nothing
end
check_plan_input(p::SeriesFFTPlan, s::AbstractSeries) = check_plan_series(input_type(p), s, p.input_metadata, "input")
check_plan_input(::SeriesFFTPlan, ::AbstractMatrix) = nothing
check_plan_destination(p::SeriesFFTPlan, s::AbstractSeries) = check_plan_series(output_type(p), s, metadata(p), "output")
check_plan_destination(::SeriesFFTPlan, ::AbstractMatrix) = nothing
storage(s::AbstractSeries) = parent(s)
storage(a::AbstractMatrix) = a

"""Apply the plan; series arguments are checked against the grids it was prepared for."""
function transform!(destination::AbstractMatrix, p::SeriesFFTPlan, input::AbstractMatrix)
    check_plan_input(p, input)
    check_plan_destination(p, destination)
    transform_storage!(storage(destination), p, storage(input))
    return destination
end
function transform(p::SeriesFFTPlan, s::AbstractSeries)
    check_plan_input(p, s)
    return output_type(p)(transform(p, parent(s)), metadata(p), Validated())
end
