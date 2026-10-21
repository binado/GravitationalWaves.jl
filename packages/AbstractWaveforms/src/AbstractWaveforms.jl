"""
    AbstractWaveforms

The interface shared by the waveform backends (`FrequencyDomainWaveforms`, `LALWaveforms`): the supertype
[`AbstractWaveformGenerator`](@ref), the generic functions [`hphc_fd`](@ref),
[`amplitude_phase_fd`](@ref) and [`validate_parameters`](@ref), and the accessors
and capability queries every generator provides. Uses the MIT GWSeries types, so generic code can depend
on it alone.
"""
module AbstractWaveforms

using GWSeries: TimeSeries, FrequencySeries, times, channels, frequencies, frequency_spacing, epoch, sampling_interval

export AbstractWaveformGenerator, AlignedSpinQuadrupolarWaveformGenerator
export HigherModesWaveformGenerator, waveform_modes_fd, waveform_modes_td, hphc_td
export hphc_fd, amplitude_phase_fd, validate_parameters
export approximant, reference_frequency, minimum_frequency, maximum_frequency
export spin_support, supports_mode_array, is_differentiable

# Helpers for backend implementers. The shared contract tests live in `test/contract.jl`, which
# backends `include` from their own tests.
public check_frequency_grid, frequency_bounds, mode_angles, mode_sum, sylm_m2

"""
    AbstractWaveformGenerator

Supertype of waveform generators in either domain. A generator holds everything about a
waveform call that does not vary between parameter draws (approximant, frequency bounds,
post-Newtonian orders, ...); build it once and reuse it with [`hphc_fd`](@ref).

A concrete subtype implements [`hphc_fd`](@ref), optionally [`validate_parameters`](@ref), and the accessors [`approximant`](@ref),
[`reference_frequency`](@ref), [`minimum_frequency`](@ref) and [`maximum_frequency`](@ref). It
may override the capability queries [`spin_support`](@ref), [`supports_mode_array`](@ref) and
[`is_differentiable`](@ref). Aligned-spin quadrupole backends can instead subtype
[`AlignedSpinQuadrupolarWaveformGenerator`](@ref) and implement [`amplitude_phase_fd`](@ref)
to inherit polarization assembly. Higher-mode backends can subtype
[`HigherModesWaveformGenerator`](@ref) and implement modes in either domain.
"""
abstract type AbstractWaveformGenerator end

"""
    AlignedSpinQuadrupolarWaveformGenerator <: AbstractWaveformGenerator

Generator for aligned-spin, quadrupole-only radiation. Implement [`amplitude_phase_fd`](@ref)
with the usual generator accessors to obtain the default [`hphc_fd`](@ref) assembly.
The inclination is `params.theta_jn`, since the orbital and total angular momenta are aligned.
Generators retain fixed configuration; parameter-dependent coefficients are evaluated per call.
"""
abstract type AlignedSpinQuadrupolarWaveformGenerator <: AbstractWaveformGenerator end

"""
    amplitude_phase_fd(f::AbstractRange{Float64}, params::NamedTuple,
                    generator::AlignedSpinQuadrupolarWaveformGenerator) -> (; amplitude, phase)

Inclination-independent frequency-domain strain represented as `amplitude .* exp.(im .* phase)`.
Both returned real arrays have `length(f)` entries. Frequencies are in Hz and phase is in radians.
Amplitude is nonnegative and includes luminosity-distance normalization. Phase is analytic and
unwrapped (not restricted to a principal interval), including orbital reference-phase alignment
and any model time shift. It is defined modulo a constant integer multiple of `2π`.

Both arrays are zero at DC, outside the configured frequency bounds, and outside intrinsic model
support. Phase is also zero wherever amplitude vanishes. This is a finite placeholder: phase
and its derivatives there have no physical meaning. Never evaluate singular formulas at DC.

Parameters, grid requirements and host validation follow [`hphc_fd`](@ref). The arrays exclude
inclination factors; changing `theta_jn` does not change their values. Differentiable backends
must support this method under the same conditions as their `hphc_fd` implementation.
"""
function amplitude_phase_fd end

"""
    hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, generator::AbstractWaveformGenerator) -> Matrix

Frequency-domain plus and cross polarisations on the uniform grid `f` (Hz).

`params` is a `NamedTuple` with bilby's parameter names and units: required `mass_1`, `mass_2`
(detector-frame solar masses, or other mass parameters that convert to them),
`luminosity_distance` (Mpc), `theta_jn`, `phase`; optional (default 0) `a_1`, `a_2`, `tilt_1`,
`tilt_2`, `phi_12`, `phi_jl`, `eccentricity`.

Whether `hphc_fd` checks the parameter values is up to the backend. Backends that can be traced
by compilers (Reactant) trust them and leave the checks to [`validate_parameters`](@ref), which
callers run once on the host. The shared polarization assemblers assume valid inputs and
backend outputs without performing runtime validation.

Grid requirements depend on the backend. The zero-start uniform-grid convention is checked
by [`check_frequency_grid`](@ref); backends may also support ranges starting above zero.
The returned matrix has `length(f)` rows and plus/cross columns and is zero outside
`minimum_frequency(generator) ≤ f ≤ maximum_frequency(generator)`, where a `nothing` maximum
means the last grid point ([`frequency_bounds`](@ref)).

For [`AlignedSpinQuadrupolarWaveformGenerator`](@ref), the default calls [`amplitude_phase_fd`](@ref)
once and multiplies `A * exp(im * Φ)` by `(1 + cos(theta_jn)^2)/2` for plus and
`-im * cos(theta_jn)` for cross.
"""
function hphc_fd end

"""Generate peak-relative polarizations on a prototype's frequency grid; its epoch and channels are ignored."""
function hphc_fd(prototype::FrequencySeries, params, g::AbstractWaveformGenerator)
    data = hphc_fd(frequencies(prototype), params, g)
    return FrequencySeries(data; epoch = 0.0, frequency_spacing = frequency_spacing(prototype), channels = (:plus, :cross))
end

function hphc_fd(f::AbstractRange{Float64}, params::NamedTuple, g::AlignedSpinQuadrupolarWaveformGenerator)
    a = amplitude_phase_fd(f, params, g)
    cfac = cos(params.theta_jn)
    pfac = 0.5 * (1.0 + cfac * cfac)
    re = a.amplitude .* cos.(a.phase)
    im = a.amplitude .* sin.(a.phase)
    plus = complex.(pfac .* re, pfac .* im)
    cross = complex.(cfac .* im, .-(cfac .* re))
    return hcat(plus, cross)
end

"""
    validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, generator::AbstractWaveformGenerator) -> nothing

Check that `hphc_fd(f, params, generator)` is well defined: the grid, the parameter names, and the
parameter values (positive masses and distance, aligned spins, ...). Throws `ArgumentError` (or
a more specific error) on invalid input and returns `nothing` otherwise.

Meant to run once on the host, before compiling `hphc_fd` or sampling with it, for the backends
whose `hphc_fd` does not check parameter values. The fallback accepts everything.
"""
validate_parameters(f::AbstractRange{Float64}, params::NamedTuple, ::AbstractWaveformGenerator) = nothing

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

Whether the generator exposes individual harmonic modes or can restrict the modes it sums.
Higher-mode generators return `true`; aligned-spin quadrupole generators return `false`.
Falls back to `false`.
"""
supports_mode_array(::AbstractWaveformGenerator) = false
supports_mode_array(::AlignedSpinQuadrupolarWaveformGenerator) = false

"""
    is_differentiable(generator) -> Bool

Whether waveform evaluation in the supported domain can be differentiated with respect to
its parameters by automatic
differentiation. Falls back to `false`.
"""
is_differentiable(::AbstractWaveformGenerator) = false

"""
    check_frequency_grid(f)

Throw `ArgumentError` unless `f` starts at 0, has a positive step and at least two points, as
bilby assumes.
"""
function check_frequency_grid(f::AbstractRange{Float64})
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
function frequency_bounds(g::AbstractWaveformGenerator, f::AbstractRange{Float64})
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

include("spherical_harmonics.jl")
include("modes.jl")

end
