"""
    TabulatedPSD(frequencies, values; kind=:psd, interpolation=:linear)

A one-sided strain PSD table, with frequencies in Hz and densities in strain²/Hz.
Use `kind=:asd` for strain/√Hz input; it is squared before interpolation.
`:linear` interpolates PSD against frequency; `:loglog` interpolates log PSD
against log frequency and requires positive table entries. Inputs are copied
to `Float64` vectors. The frequency table may be nonuniform.
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

"""
    read_psd(path_or_io; interpolation=:linear)

Read a two-column frequency/PSD table. Whitespace and comma delimiters, blank
lines and `#` comments are supported. No header or unit conversion is inferred.
"""
read_psd(source; interpolation = :linear) = read_density(source, :psd, interpolation)

"""
    read_asd(path_or_io; interpolation=:linear)

Read a two-column frequency/ASD table and convert to a `TabulatedPSD`.
The text format is the same as [`read_psd`](@ref).
"""
read_asd(source; interpolation = :linear) = read_density(source, :asd, interpolation)

function read_density(path::AbstractString, kind, interpolation)
    return open(io -> read_density(io, kind, interpolation), path)
end

function read_density(io::IO, kind, interpolation)
    text = read(io, String)
    row = first_data_row(text)
    isnothing(row) && throw(ArgumentError("no frequency/density rows found"))
    # The delimiter is decided by the first data row, since comments may contain commas.
    table = if occursin(',', row)
        readdlm(IOBuffer(text), ',', Float64; comments = true)
    else
        readdlm(IOBuffer(text), Float64; comments = true)
    end::Matrix{Float64}
    size(table, 2) == 2 || throw(ArgumentError("expected two columns: frequency and density"))
    return TabulatedPSD(table[:, 1], table[:, 2]; kind, interpolation)
end

function first_data_row(text)
    for line in eachline(IOBuffer(text))
        row = strip(first(split(line, '#'; limit = 2)))
        isempty(row) || return row
    end
    return nothing
end

"""
    psd(f, curve::TabulatedPSD; outside=:error)

Evaluate a one-sided PSD at a frequency or collection of frequencies (Hz).
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

psd(f, curve::TabulatedPSD; outside = :error) = map(x -> psd(x, curve; outside), f)

"""`asd(f, curve; outside=:error)` evaluates the square root of [`psd`](@ref)."""
asd(f, curve::TabulatedPSD; outside = :error) = sqrt.(psd(f, curve; outside))
