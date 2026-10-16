"""
    Approximant

Supertype of the waveform models implemented here. `WaveformSettings{A}` is parametric in the
model so `hphc_fd` dispatches on it.
"""
abstract type Approximant end

"Post-Newtonian stationary-phase inspiral (`XLALSimInspiralTaylorF2`), aligned spins only."
struct TaylorF2 <: Approximant end

const APPROXIMANTS = (TaylorF2 = TaylorF2(),)

function approximant_from_name(name::Union{Symbol, AbstractString})
    key = Symbol(name)
    haskey(APPROXIMANTS, key) || throw(
        ArgumentError(
            "unknown approximant `$name`; implemented approximants are " *
                join(keys(APPROXIMANTS), ", ")
        )
    )
    return APPROXIMANTS[key]
end

approximant_name(::A) where {A <: Approximant} = String(nameof(A))
