# Experimental only: all precomputed physical coefficients remain tunable.
module StructuredParametersPrototype
using GWSEOBNR
import SciMLStructures as SS
const S = GWSEOBNR
const C_NAMES = (:X_1, :X_2, :chi_1, :chi_2, :nu, :a6, :dSO)
const F_NAMES = (:rho, :rholog, :f, :fvh, :delta, :deltavh, :extra, :extralog)
const BLOCK = 9 * 9 * 11
const COUNT = 7 + 8BLOCK + 2 * 81

struct EOBParameters{C, F}
    c::C
    coeffs::F
end
SS.isscimlstructure(::EOBParameters) = true
SS.ismutablescimlstructure(::EOBParameters) = false
SS.hasportion(::SS.AbstractPortion, ::EOBParameters) = false
SS.hasportion(::SS.Tunable, ::EOBParameters) = true
SS.canonicalize(::SS.AbstractPortion, ::EOBParameters) = (nothing, nothing, nothing)

function SS.canonicalize(::SS.Tunable, p::EOBParameters)
    values = vcat(
        collect(values_of_c(p.c)), vec(p.coeffs.rho), vec(p.coeffs.rholog), vec(p.coeffs.f), vec(p.coeffs.fvh),
        vec(p.coeffs.delta), vec(p.coeffs.deltavh), vec(p.coeffs.extra), vec(p.coeffs.extralog),
        vec(real.(p.coeffs.prefixes)), vec(imag.(p.coeffs.prefixes))
    )
    return values, x -> SS.replace(SS.Tunable(), p, x), false
end
values_of_c(c) = values(c)
function SS.replace(::SS.Tunable, ::EOBParameters, x)
    length(x) == COUNT || throw(DimensionMismatch("incorrect EOB tunable count"))
    c = NamedTuple{C_NAMES}(ntuple(i -> x[i], 7))
    blocks = ntuple(i -> reshape(x[(8 + (i - 1) * BLOCK):(7 + i * BLOCK)], 9, 9, 11), 8)
    offset = 7 + 8BLOCK
    prefixes = reshape(complex.(x[(offset + 1):(offset + 81)], x[(offset + 82):(offset + 162)]), 9, 9)
    coeffs = merge(NamedTuple{F_NAMES}(blocks), (; prefixes))
    return EOBParameters(c, coeffs)
end

# Install only in a dedicated exploratory process. Physics/solver/timing code is
# loaded verbatim from production; the sole substitution is the ODE container.
function install!()
    path = joinpath(pkgdir(S), "src", "dynamics.jl")
    source = read(path, String)
    source = source[findfirst("function evolve(", source).start:end]
    needle = "pars = (; c, coeffs)"
    occursin(needle, source) || error("production parameter construction changed")
    Core.eval(S, :(const ExperimentalEOBParameters = $EOBParameters))
    Base.include_string(S, replace(source, needle => "pars = ExperimentalEOBParameters(c, coeffs)"), path)
    return nothing
end
end
