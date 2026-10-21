function check_fourier_domain(s::TimeSeries)
    nsamples(s) >= 2 || throw(ArgumentError("FFT needs at least two time samples"))
    return nothing
end
function check_fourier_domain(s::FrequencySeries; nsamples = nothing)
    if metadata(s).spectrum == :onesided
        n = something(nsamples, 2 * (size(s, 1) - 1))
        n >= 2 && size(s, 1) == n ÷ 2 + 1 && first(frequencies(s)) == 0 ||
            throw(ArgumentError("inverse real FFT needs the complete zero-origin one-sided grid for its output length"))
    else
        n = size(s, 1)
        n >= 2 && (nsamples === nothing || nsamples == n) || throw(ArgumentError("inverse FFT length differs"))
        check_coordinate(first(frequencies(s)), -fld(n, 2) * frequency_spacing(s), frequency_spacing(s), "full-spectrum origins")
    end
    return nothing
end

function transform_dims(dims)
    (dims == 1 || dims == (1,) || dims == [1]) || throw(ArgumentError("series FFTs transform only the sample dimension (1)"))
    return 1
end

function fourier_metadata(s::TimeSeries, real_transform)
    df = inv(duration(s))
    return (;
        epoch = epoch(s), frequency_start = real_transform ? zero(df) : -fld(nsamples(s), 2) * df,
        frequency_spacing = df, channels = channels(s), spectrum = real_transform ? :onesided : :twosided,
    )
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
    return FrequencySeries(a, fourier_metadata(s, true), Validated())
end
"""Physical dt*fft along rows, stored with increasing centered signed frequencies."""
function AbstractFFTs.fft(s::TimeSeries, dims = 1)
    transform_dims(dims)
    check_fourier_domain(s)
    a = AbstractFFTs.fftshift(AbstractFFTs.fft(parent(s), 1), 1)
    a .*= sampling_interval(s)
    return FrequencySeries(a, fourier_metadata(s, false), Validated())
end
function AbstractFFTs.irfft(s::FrequencySeries, n::Integer, dims = 1)
    transform_dims(dims)
    metadata(s).spectrum == :onesided || throw(ArgumentError("irfft requires a one-sided spectrum"))
    check_fourier_domain(s; nsamples = n)
    a = AbstractFFTs.irfft(complex.(parent(s)), n, 1)
    a .*= n * frequency_spacing(s)
    return TimeSeries(a, inverse_metadata(s, n), Validated())
end
function AbstractFFTs.ifft(s::FrequencySeries, dims = 1)
    transform_dims(dims)
    metadata(s).spectrum == :twosided || throw(ArgumentError("ifft requires a centered two-sided spectrum; use irfft for one-sided data"))
    check_fourier_domain(s)
    a = AbstractFFTs.ifft(AbstractFFTs.ifftshift(parent(s), 1), 1)
    a .*= nsamples(s) * frequency_spacing(s)
    return TimeSeries(a, inverse_metadata(s, nsamples(s)), Validated())
end
for f in (:fft!, :ifft!, :bfft!, :plan_fft!, :plan_ifft!, :plan_bfft!)
    @eval AbstractFFTs.$f(s::AbstractSeries, args...; kwargs...) = throw(ArgumentError("in-place FFTs cannot change a series type; use a physical transform plan or parent(series)"))
end
for f in (:fftshift, :ifftshift)
    @eval AbstractFFTs.$f(s::AbstractSeries, args...) = throw(ArgumentError("shifting bins changes series coordinates; use parent(series) or the physical FFT methods"))
end

"""Transform a time series with physical normalization (real FFT for real samples)."""
to_frequencydomain(s::TimeSeries{<:Real}) = AbstractFFTs.rfft(s)
to_frequencydomain(s::TimeSeries) = AbstractFFTs.fft(s)
"""Invert a physical spectrum; specify nsamples for an odd real-transform length."""
function to_timedomain(s::FrequencySeries; nsamples = nothing)
    if metadata(s).spectrum == :onesided
        return AbstractFFTs.irfft(s, something(nsamples, 2 * (size(s, 1) - 1)))
    end
    check_fourier_domain(s; nsamples)
    return AbstractFFTs.ifft(s)
end

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
function plan_to_frequencydomain(s::TimeSeries; real_transform = eltype(s) <: Real, kwargs...)
    check_fourier_domain(s)
    real_transform && !(eltype(s) <: Real) && throw(ArgumentError("real FFT needs real samples"))
    input = real_transform ? copy(parent(s)) : complex.(parent(s))
    fft = real_transform ? AbstractFFTs.plan_rfft(input, 1; kwargs...) : AbstractFFTs.plan_fft(input, 1; kwargs...)
    shape = (real_transform ? nsamples(s) ÷ 2 + 1 : nsamples(s), size(s, 2))
    work = similar(input, Complex{real(eltype(input))}, shape)
    return SeriesFFTPlan(fft, input, work, metadata(s), fourier_metadata(s, real_transform), sampling_interval(s), real_transform ? :rfft : :fft)
end
"""Prepare a source-preserving physical inverse FFT, optionally specifying real output length."""
function plan_to_timedomain(s::FrequencySeries; nsamples = nothing, kwargs...)
    check_fourier_domain(s; nsamples)
    real_transform = metadata(s).spectrum == :onesided
    n = real_transform ? something(nsamples, 2 * (size(s, 1) - 1)) : size(s, 1)
    input = complex.(parent(s))
    fft = real_transform ? AbstractFFTs.plan_irfft(input, n, 1; kwargs...) : AbstractFFTs.plan_ifft(input, 1; kwargs...)
    T = real_transform ? real(eltype(input)) : eltype(input)
    work = similar(input, T, (n, size(s, 2)))
    return SeriesFFTPlan(fft, input, work, metadata(s), inverse_metadata(s, n), n * frequency_spacing(s), real_transform ? :irfft : :ifft)
end

function transform!(destination::AbstractMatrix, p::SeriesFFTPlan, input::AbstractMatrix)
    Base.mightalias(destination, input) && throw(ArgumentError("transform destination must not alias its input"))
    size(input) == size(p.input) && size(destination) == size(p.work) || throw(DimensionMismatch("transform storage does not match the prepared shapes"))
    if p.kind == :ifft
        AbstractFFTs.ifftshift!(p.input, input, 1)
    else
        copyto!(p.input, input)
    end
    mul!(p.work, p.fft, p.input)
    if p.kind == :fft
        AbstractFFTs.fftshift!(destination, p.work, 1)
    else
        copyto!(destination, p.work)
    end
    destination .*= p.scale
    return destination
end
function transform(p::SeriesFFTPlan, input::AbstractMatrix)
    return transform!(similar(p.work), p, input)
end
function checked_plan_input(p, s)
    forward = p.kind in (:fft, :rfft)
    (forward ? s isa TimeSeries : s isa FrequencySeries) || throw(ArgumentError("series kind does not match transform input"))
    prototype = forward ? TimeSeries(parent(s); p.input_metadata...) : FrequencySeries(parent(s); p.input_metadata...)
    check_same_grid(prototype, s)
    return nothing
end
function checked_plan_destination(p, s)
    forward = p.kind in (:fft, :rfft)
    (forward ? s isa FrequencySeries : s isa TimeSeries) || throw(ArgumentError("series kind does not match transform output"))
    expected = forward ? FrequencySeries(parent(s); metadata(p)...) : TimeSeries(parent(s); metadata(p)...)
    check_same_grid(expected, s)
    return nothing
end
function transform(p::SeriesFFTPlan, s::AbstractSeries)
    checked_plan_input(p, s)
    a = transform(p, parent(s))
    return p.kind in (:fft, :rfft) ? FrequencySeries(a; metadata(p)...) : TimeSeries(a; metadata(p)...)
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
