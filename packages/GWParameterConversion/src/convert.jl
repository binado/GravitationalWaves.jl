# Key-driven conversion between parameter NamedTuples.
#
# Every conversion is a registry entry `needs -> gives` with a plain function. For a given set
# of input keys and requested targets the plan (which entries to run, in which order) is
# resolved at compile time by a generated function, so a call compiles down to the needed
# arithmetic and is type stable.

struct Conversion{N, M, F}
    needs::NTuple{N, Symbol}
    gives::NTuple{M, Symbol}
    f::F
end

# Order matters: when several routes could fill a key, the earliest registered entry wins, and
# component masses are always preferred as the source for the other mass parameters.
const MASS_CONVERSIONS = (
    Conversion((:mass_1, :mass_2), (:chirp_mass,), chirp_mass),
    Conversion((:mass_1, :mass_2), (:total_mass,), total_mass),
    Conversion((:mass_1, :mass_2), (:mass_ratio,), mass_ratio),
    Conversion((:mass_1, :mass_2), (:symmetric_mass_ratio,), symmetric_mass_ratio),
    Conversion(
        (:chirp_mass, :mass_ratio), (:mass_1, :mass_2),
        chirp_mass_and_mass_ratio_to_component_masses,
    ),
    Conversion(
        (:total_mass, :mass_ratio), (:mass_1, :mass_2),
        total_mass_and_mass_ratio_to_component_masses,
    ),
    Conversion((:symmetric_mass_ratio,), (:mass_ratio,), symmetric_mass_ratio_to_mass_ratio),
    Conversion(
        (:chirp_mass, :total_mass), (:symmetric_mass_ratio,),
        chirp_mass_and_total_mass_to_symmetric_mass_ratio,
    ),
    Conversion(
        (:chirp_mass, :mass_ratio), (:total_mass,), chirp_mass_and_mass_ratio_to_total_mass
    ),
    Conversion((:mass_1, :chirp_mass), (:mass_ratio,), mass_1_and_chirp_mass_to_mass_ratio),
    Conversion((:mass_2, :chirp_mass), (:mass_ratio,), mass_2_and_chirp_mass_to_mass_ratio),
    # Source <-> detector frame (the redshift is an input; deriving it needs a cosmology).
    (
        Conversion((name, :redshift), (Symbol(name, :_source),), detector_to_source_mass)
            for name in (:mass_1, :mass_2, :chirp_mass, :total_mass)
    )...,
    (
        Conversion((Symbol(name, :_source), :redshift), (name,), source_to_detector_mass)
            for name in (:mass_1, :mass_2, :chirp_mass, :total_mass)
    )...,
)

const SPIN_CONVERSIONS = (
    Conversion((:a_1, :tilt_1), (:spin_1z,), (a, tilt) -> a * cos(tilt)),
    Conversion((:a_2, :tilt_2), (:spin_2z,), (a, tilt) -> a * cos(tilt)),
    Conversion((:a_1, :tilt_1), (:chi_1_in_plane,), (a, tilt) -> a * sin(tilt)),
    Conversion((:a_2, :tilt_2), (:chi_2_in_plane,), (a, tilt) -> a * sin(tilt)),
    Conversion((:spin_1z, :spin_2z, :mass_ratio), (:chi_eff,), effective_spin),
    Conversion(
        (:chi_1_in_plane, :chi_2_in_plane, :mass_ratio), (:chi_p,), effective_precessing_spin
    ),
)

const REGISTRY = (MASS_CONVERSIONS..., SPIN_CONVERSIONS...)

"""
    MASS_KEYS

Every parameter name that takes part in a mass conversion.
"""
const MASS_KEYS = Tuple(
    unique(Iterators.flatten((Iterators.flatten((c.needs, c.gives)) for c in MASS_CONVERSIONS)))
)

"""
    InconsistentParametersError(key, given, derived, from)

Thrown by [`convert_parameters`](@ref) when `key` was supplied together with the parameters
`from` that determine it, and the two disagree.
"""
struct InconsistentParametersError <: Exception
    key::Symbol
    given::Any
    derived::Any
    from::Tuple{Vararg{Symbol}}
end

function Base.showerror(io::IO, e::InconsistentParametersError)
    return print(
        io, "InconsistentParametersError: `", e.key, " = ", e.given, "` disagrees with the value ",
        e.derived, " implied by ", join(e.from, ", ")
    )
end

_tuple(x::Tuple) = x
_tuple(x) = (x,)

function _check(given, derived, key, from, rtol)
    isapprox(given, derived; rtol, atol = rtol) ||
        throw(InconsistentParametersError(key, given, derived, from))
    return nothing
end

# Plan the entries to run: a forward pass applies every entry that can fill something new until
# the targets are available (or nothing changes), then a backward pass keeps only the entries
# needed for the targets. Returns the entry indices, or an error message.
function _plan(keys::Tuple, targets::Tuple)
    available = Set{Symbol}(keys)
    applied = Int[]
    progress = true
    while progress && !all(in(available), targets)
        progress = false
        for (i, c) in enumerate(REGISTRY)
            all(in(available), c.needs) || continue
            all(in(available), c.gives) && continue
            push!(applied, i)
            union!(available, c.gives)
            progress = true
        end
    end
    unreachable = filter(!in(available), collect(targets))
    isempty(unreachable) ||
        return "cannot derive $(join(unreachable, ", ")) from the keys ($(join(keys, ", ")))"

    given = Set{Symbol}(keys)
    needed = Set{Symbol}(filter(!in(given), collect(targets)))
    keep = Int[]
    for i in reverse(applied)
        c = REGISTRY[i]
        any(in(needed), c.gives) || continue
        pushfirst!(keep, i)
        union!(needed, filter(!in(given), collect(c.needs)))
    end
    return keep
end

@generated function _convert(nt::NamedTuple{K}, ::Val{T}, ::Val{C}, rtol) where {K, T, C}
    plan = _plan(K, T)
    plan isa String && return :(throw(ArgumentError($plan)))

    body = Expr[]
    # Redundant inputs: check them against what the other given keys imply.
    for c in (C ? REGISTRY : ())
        all(in(K), c.needs) || continue
        present = [j for (j, g) in enumerate(c.gives) if g in K]
        isempty(present) && continue
        args = [:(getfield(nt, $(QuoteNode(n)))) for n in c.needs]
        checks = [
            :(
                _check(
                    getfield(nt, $(QuoteNode(c.gives[j]))), vals[$j], $(QuoteNode(c.gives[j])),
                    $(QuoteNode(c.needs)), rtol
                )
            )
                for j in present
        ]
        push!(
            body, :(
                let vals = _tuple($(c.f)($(args...)))
                    $(checks...)
                end
            )
        )
    end

    available = Set{Symbol}(K)
    push!(body, :(acc = nt))
    for i in plan
        c = REGISTRY[i]
        js = [j for (j, g) in enumerate(c.gives) if !(g in available)]
        isempty(js) && continue
        args = [:(getfield(acc, $(QuoteNode(n)))) for n in c.needs]
        names = Expr(:tuple, (QuoteNode(c.gives[j]) for j in js)...)
        values = Expr(:tuple, (:(vals[$j]) for j in js)...)
        push!(
            body, :(
                let vals = _tuple($(c.f)($(args...)))
                    acc = merge(acc, NamedTuple{$names}($values))
                end
            )
        )
        union!(available, c.gives)
    end
    return Expr(:block, body..., :(acc))
end

"""
    convert_parameters(nt::NamedTuple, targets::Tuple{Vararg{Symbol}}; rtol = 1e-12, check_consistency::Bool = true) -> NamedTuple

Return `nt` extended with the parameters in `targets` that it lacks, derived from the keys it
has (for example `(:mass_1, :mass_2)` from `chirp_mass` and `mass_ratio`). Keys already in `nt`
are never overwritten, and unrelated keys pass through.

Conversions covered: all mass parameters (`mass_1`, `mass_2`, `chirp_mass`, `total_mass`,
`mass_ratio`, `symmetric_mass_ratio`), their `_source` versions given a `redshift`, and
`spin_1z`, `spin_2z`, `chi_1_in_plane`, `chi_2_in_plane`, `chi_eff`, `chi_p` from magnitudes and
tilts. The route is resolved at compile time from the key names, provided `targets` is a literal
or a constant (otherwise the call still works but is not type stable); when several routes exist
the earliest registered one is used (component masses are the preferred source).

Throws `ArgumentError` when a target cannot be derived, and [`InconsistentParametersError`](@ref)
when `nt` holds redundant parameters that disagree, i.e. differ by more than `rtol` (relative and
absolute). Set `check_consistency=false` for compiled inference calls: redundant-input
calculations and checks are omitted entirely, so callers must supply consistent redundant
inputs. Supplied keys still win. Keep `targets`, `check_consistency`, and `rtol` as
compile-time configuration when tracing.
"""
@inline function convert_parameters(
        nt::NamedTuple, targets::Tuple{Vararg{Symbol}}; rtol::Number = 1.0e-12,
        check_consistency::Bool = true
    )
    return _convert(nt, Val(targets), Val(check_consistency), rtol)
end
