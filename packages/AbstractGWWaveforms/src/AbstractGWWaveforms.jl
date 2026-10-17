"""
    AbstractGWWaveforms

The interface shared by the waveform backends (`GWWaveforms`, `LALWaveforms`): the supertype
[`AbstractWaveformGenerator`](@ref), the generic functions [`hphc_fd`](@ref) and
[`validate_parameters`](@ref) and the accessors
and capability queries every generator provides. Dependency-free, so generic code can depend
on it alone.
"""
module AbstractGWWaveforms

export AbstractWaveformGenerator, hphc_fd, validate_parameters
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable

# Helpers for backend implementers. The shared contract tests live in `test/contract.jl`, which
# backends `include` from their own tests.
public check_frequency_grid, frequency_bounds

"""
    AbstractWaveformGenerator

Supertype of frequency-domain waveform generators. A generator holds everything about a
waveform call that does not vary between parameter draws (approximant, frequency bounds,
post-Newtonian orders, ...); build it once and reuse it with [`hphc_fd`](@ref).

A concrete subtype implements [`hphc_fd`](@ref), optionally [`validate_parameters`](@ref), and the accessors [`approximant`](@ref),
[`reference_frequency`](@ref), [`minimum_frequency`](@ref) and [`maximum_frequency`](@ref). It
may override the capability queries [`spin_support`](@ref), [`supports_mode_array`](@ref) and
[`is_differentiable`](@ref).
"""
abstract type AbstractWaveformGenerator end

"""
    hphc_fd(f::AbstractRange, params::NamedTuple, generator::AbstractWaveformGenerator) -> (; plus, cross)

Frequency-domain plus and cross polarisations on the uniform grid `f` (Hz).

`params` is a `NamedTuple` with bilby's parameter names and units: required `mass_1`, `mass_2`
(detector-frame solar masses, or other mass parameters that convert to them),
`luminosity_distance` (Mpc), `theta_jn`, `phase`; optional (default 0) `a_1`, `a_2`, `tilt_1`,
`tilt_2`, `phi_12`, `phi_jl`, `eccentricity`.

Whether `hphc_fd` checks the parameter values is up to the backend. Backends that can be traced
by compilers (Reactant) trust them and leave the checks to [`validate_parameters`](@ref), which
callers run once on the host. Properties of the grid and the generator stay checked.

`f` must start at 0, have a positive step and at least two points ([`check_frequency_grid`](@ref)).
The returned arrays have `length(f)` entries and are zero outside
`minimum_frequency(generator) ≤ f ≤ maximum_frequency(generator)`, where a `nothing` maximum
means the last grid point ([`frequency_bounds`](@ref)).
"""
function hphc_fd end

"""
    validate_parameters(f::AbstractRange, params::NamedTuple, generator::AbstractWaveformGenerator) -> nothing

Check that `hphc_fd(f, params, generator)` is well defined: the grid, the parameter names, and the
parameter values (positive masses and distance, aligned spins, ...). Throws `ArgumentError` (or
a more specific error) on invalid input and returns `nothing` otherwise.

Meant to run once on the host, before compiling `hphc_fd` or sampling with it, for the backends
whose `hphc_fd` does not check parameter values. The fallback accepts everything.
"""
validate_parameters(f::AbstractRange, params::NamedTuple, ::AbstractWaveformGenerator) = nothing

"""
    approximant(generator) -> String

Name of the waveform model, as LAL spells it.
"""
function approximant end

"""
    reference_frequency(generator)

Frequency (Hz) at which the orbital phase and spin angles are defined.
"""
function reference_frequency end

"""
    minimum_frequency(generator)

Lower edge (Hz) of the support of the returned waveform.
"""
function minimum_frequency end

"""
    maximum_frequency(generator)

Upper edge (Hz) of the support of the returned waveform, or `nothing` for the last point of
the frequency grid.
"""
function maximum_frequency end

"""
    spin_support(generator) -> Symbol

Spin degrees of freedom the model evolves: `:none`, `:aligned` (spins along the orbital angular
momentum) or `:precessing`. Falls back to `:aligned`.
"""
spin_support(::AbstractWaveformGenerator) = :aligned

"""
    supports_mode_array(generator) -> Bool

Whether the generator can restrict the harmonic modes it sums. Falls back to `false`.
"""
supports_mode_array(::AbstractWaveformGenerator) = false

"""
    is_differentiable(generator) -> Bool

Whether `hphc_fd` can be differentiated with respect to its parameters by automatic
differentiation. Falls back to `false`.
"""
is_differentiable(::AbstractWaveformGenerator) = false

"""
    check_frequency_grid(f)

Throw `ArgumentError` unless `f` starts at 0, has a positive step and at least two points, as
bilby assumes.
"""
function check_frequency_grid(f::AbstractRange)
    first(f) == 0 || throw(ArgumentError("frequency grid must start at 0, got $(first(f))"))
    df = step(f)
    df > 0 || throw(ArgumentError("frequency grid must have a positive step, got $df"))
    length(f) >= 2 || throw(ArgumentError("frequency grid needs at least two points"))
    return f
end

"""
    frequency_bounds(generator, f) -> (f_min, f_max)

Support `[f_min, f_max]` of the waveform on the grid `f`.
"""
function frequency_bounds(g::AbstractWaveformGenerator, f::AbstractRange)
    return (minimum_frequency(g), something(maximum_frequency(g), Float64(last(f))))
end

function Base.show(io::IO, g::AbstractWaveformGenerator)
    return print(
        io, nameof(typeof(g)), "(", repr(approximant(g)),
        "; reference_frequency=", reference_frequency(g),
        ", minimum_frequency=", minimum_frequency(g),
        ", maximum_frequency=", maximum_frequency(g), ")"
    )
end

end
