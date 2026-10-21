function check_fourier_domain(s::TimeSeries)
    nsamples(s) >= 2 || throw(ArgumentError("FFT needs at least two time samples"))
    return nothing
end
function check_fourier_domain(s::FrequencySeries; nsamples::Integer)
    nsamples >= 2 && size(s, 1) == nsamples ÷ 2 + 1 ||
        throw(ArgumentError("inverse real FFT needs N÷2+1 bins for its output length"))
    return nothing
end

function transform_dims(dims)
    (dims == 1 || dims == (1,) || dims == [1]) || throw(ArgumentError("series FFTs transform only the sample dimension (1)"))
    return 1
end

function fourier_metadata(s::TimeSeries)
    return (; epoch = epoch(s), frequency_spacing = inv(duration(s)), channels = channels(s))
end
function inverse_metadata(s::FrequencySeries, n)
    return (; epoch = epoch(s), sampling_interval = inv(n * frequency_spacing(s)), channels = channels(s))
end

"""Physical dt*rfft along rows; the loaded AbstractFFTs backend owns execution."""
function AbstractFFTs.rfft(s::TimeSeries{<:Real}, dims = 1)
    transform_dims(dims)
    check_fourier_domain(s)
    a = AbstractFFTs.rfft(parent(s), 1)
    a .*= sampling_interval(s)
    return FrequencySeries(a, fourier_metadata(s), Validated())
end
AbstractFFTs.rfft(s::TimeSeries, dims = 1) = throw(ArgumentError("wrapped forward transforms require real samples"))

function AbstractFFTs.irfft(s::FrequencySeries, n::Integer, dims = 1)
    transform_dims(dims)
    check_fourier_domain(s; nsamples = n)
    a = AbstractFFTs.irfft(complex.(parent(s)), n, 1)
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
"""Invert a physical one-sided spectrum with explicit original time length."""
to_timedomain(s::FrequencySeries; nsamples::Integer) = AbstractFFTs.irfft(s, nsamples)

"""Reusable AbstractFFTs plan and private scratch storage; not reentrant."""
struct SeriesFFTPlan{P, I, W, IM, OM, S}
    fft::P
    input::I
    work::W
    input_metadata::IM
    output_metadata::OM
    scale::S
    kind::Symbol
end
metadata(p::SeriesFFTPlan) = p.output_metadata

"""Prepare a physical forward FFT. Backend-specific keywords pass to AbstractFFTs."""
function plan_to_frequencydomain(s::TimeSeries{<:Real}; kwargs...)
    check_fourier_domain(s)
    input = copy(parent(s))
    fft = AbstractFFTs.plan_rfft(input, 1; kwargs...)
    shape = (nsamples(s) ÷ 2 + 1, size(s, 2))
    work = similar(input, Complex{real(eltype(input))}, shape)
    return SeriesFFTPlan(fft, input, work, metadata(s), fourier_metadata(s), sampling_interval(s), :rfft)
end
plan_to_frequencydomain(s::TimeSeries; kwargs...) = throw(ArgumentError("wrapped forward transforms require real samples"))

"""Prepare a source-preserving inverse real FFT with explicit original time length."""
function plan_to_timedomain(s::FrequencySeries; nsamples::Integer, kwargs...)
    check_fourier_domain(s; nsamples)
    input = complex.(parent(s))
    fft = AbstractFFTs.plan_irfft(input, nsamples, 1; kwargs...)
    work = similar(input, real(eltype(input)), (nsamples, size(s, 2)))
    return SeriesFFTPlan(fft, input, work, metadata(s), inverse_metadata(s, nsamples), nsamples * frequency_spacing(s), :irfft)
end

function transform!(destination::AbstractMatrix, p::SeriesFFTPlan, input::AbstractMatrix)
    Base.mightalias(destination, input) && throw(ArgumentError("transform destination must not alias its input"))
    size(input) == size(p.input) && size(destination) == size(p.work) || throw(DimensionMismatch("transform storage does not match the prepared shapes"))
    copyto!(p.input, input)
    mul!(p.work, p.fft, p.input)
    destination .= p.work .* p.scale
    return destination
end
function transform(p::SeriesFFTPlan, input::AbstractMatrix)
    return transform!(similar(p.work), p, input)
end
function checked_plan_input(p, s)
    forward = p.kind == :rfft
    (forward ? s isa TimeSeries : s isa FrequencySeries) || throw(ArgumentError("series kind does not match transform input"))
    prototype = forward ? TimeSeries(parent(s); p.input_metadata...) : FrequencySeries(parent(s); p.input_metadata...)
    check_same_grid(prototype, s)
    return nothing
end
function checked_plan_destination(p, s)
    forward = p.kind == :rfft
    (forward ? s isa FrequencySeries : s isa TimeSeries) || throw(ArgumentError("series kind does not match transform output"))
    expected = forward ? FrequencySeries(parent(s); metadata(p)...) : TimeSeries(parent(s); metadata(p)...)
    check_same_grid(expected, s)
    return nothing
end
function transform(p::SeriesFFTPlan, s::AbstractSeries)
    checked_plan_input(p, s)
    a = transform(p, parent(s))
    return p.kind == :rfft ? FrequencySeries(a; metadata(p)...) : TimeSeries(a; metadata(p)...)
end
function transform!(destination::AbstractSeries, p::SeriesFFTPlan, input::AbstractSeries)
    checked_plan_input(p, input)
    checked_plan_destination(p, destination)
    transform!(parent(destination), p, parent(input))
    return destination
end
function transform!(destination::AbstractMatrix, p::SeriesFFTPlan, input::AbstractSeries)
    checked_plan_input(p, input)
    transform!(destination, p, parent(input))
    return destination
end
function transform!(destination::AbstractSeries, p::SeriesFFTPlan, input::AbstractMatrix)
    checked_plan_destination(p, destination)
    transform!(parent(destination), p, input)
    return destination
end
