"""
    TabulatedPSD(frequencies, values; kind=:psd, interpolation=:linear)
    TabulatedPSD(table::AbstractMatrix; kind=:psd, interpolation=:linear)
    TabulatedPSD(path_or_io, delim::AbstractChar; kind=:psd, interpolation=:linear)

A one-sided strain PSD table, with frequencies in Hz and densities in strain²/Hz.
Use `kind=:asd` for strain/√Hz input; it is squared before interpolation.
`:linear` interpolates PSD against frequency; `:loglog` interpolates log PSD
against log frequency and requires positive table entries. Inputs are copied
to `Float64` vectors. The frequency table may be nonuniform.

Matrices and files must have exactly two columns: frequency and density.
File delimiters are passed directly to `DelimitedFiles.readdlm`; blank lines
and `#` comments are supported. No delimiter, header or unit conversion is inferred.
"""
struct TabulatedPSD
    frequencies::Vector{Float64}
    values::Vector{Float64}
    interpolation::Symbol

    function TabulatedPSD(frequencies, values; kind = :psd, interpolation = :linear)
        length(frequencies) == length(values) >= 2 ||
            throw(ArgumentError("need at least two matching frequencies and densities"))
        kind in (:psd, :asd) || throw(ArgumentError("kind must be :psd or :asd"))
        interpolation in (:linear, :loglog) ||
            throw(ArgumentError("interpolation must be :linear or :loglog"))
        f, s = collect(Float64, frequencies), collect(Float64, values)
        all(isfinite, f) && all(>=(0), f) && all(>(0), diff(f)) ||
            throw(ArgumentError("frequencies must be finite, nonnegative and strictly increasing"))
        all(isfinite, s) && all(>=(0), s) ||
            throw(ArgumentError("densities must be finite and nonnegative"))
        kind === :asd && (s .^= 2)
        if interpolation === :loglog
            first(f) > 0 && all(>(0), s) ||
                throw(ArgumentError("loglog interpolation requires positive frequencies and densities"))
        end
        return new(f, s, interpolation)
    end
end

Base.broadcastable(curve::TabulatedPSD) = Ref(curve)

function TabulatedPSD(table::AbstractMatrix; kind = :psd, interpolation = :linear)
    size(table, 2) == 2 || throw(ArgumentError("expected two columns: frequency and density"))
    return TabulatedPSD(table[:, 1], table[:, 2]; kind, interpolation)
end

function TabulatedPSD(
        source::Union{AbstractString, IO}, delim::AbstractChar;
        kind = :psd, interpolation = :linear
    )
    table = readdlm(source, delim, Float64; comments = true)::Matrix{Float64}
    return TabulatedPSD(table; kind, interpolation)
end

"""
    psd(f, curve::TabulatedPSD; outside=:error)

Evaluate a one-sided PSD at a frequency (Hz). Use `psd.(frequencies, curve)`
for a collection of frequencies; the curve broadcasts as a scalar.
Frequencies outside the table throw unless `outside=:zero` is explicit.
"""
function psd(f::Real, curve::TabulatedPSD; outside = :error)
    outside in (:error, :zero) || throw(ArgumentError("outside must be :error or :zero"))
    x, y = curve.frequencies, curve.values
    if !(first(x) <= f <= last(x))
        outside === :zero && return 0.0
        throw(DomainError(f, "frequency outside PSD coverage"))
    end
    i = min(searchsortedlast(x, f), length(x) - 1)
    if curve.interpolation === :loglog
        w = log(f / x[i]) / log(x[i + 1] / x[i])
        return exp((1 - w) * log(y[i]) + w * log(y[i + 1]))
    end
    w = (f - x[i]) / (x[i + 1] - x[i])
    return (1 - w) * y[i] + w * y[i + 1]
end

"""`asd(f, curve; outside=:error)` evaluates the square root of [`psd`](@ref)."""
asd(f::Real, curve::TabulatedPSD; outside = :error) = sqrt(psd(f, curve; outside))
