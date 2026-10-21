"""
    SurrogatePyCallWaveforms

Experimental PyCall adapter around upstream `gwsurrogate` implementing the
`AbstractWaveforms` time-domain higher-mode interface. Configure PyCall with the
root uv project's `gwsurrogate` extra before building PyCall.
"""
module SurrogatePyCallWaveforms

using PyCall
import AbstractWaveforms
using AbstractWaveforms: HigherModesWaveformGenerator
import AbstractWaveforms: waveform_modes_td, mode_angles

export PyCallGenerator, load_generator, waveform_modes_td

const MSUN_IN_SECONDS = 4.925491025543576e-6
const SPEED_OF_LIGHT = 299792458.0
const PARSEC = 3.085677581491367e16

"""A fixed configuration wrapper for an upstream NRSur7dq4 model."""
struct PyCallGenerator{P} <: HigherModesWaveformGenerator
    model::P
    name::String
    ell_max::Int
    f_ref::Float64
    f_low::Float64
end

"""
    load_generator(path; model=:NRSur7dq4, ell_max=4, f_ref, f_low=0)

Load an upstream `gwsurrogate` model from an existing HDF5 file. Pass the same file
path used by `SurrogateWaveforms.data_path` to share the Julia package's cached model data.
"""
function load_generator(path::AbstractString; model::Symbol = :NRSur7dq4, ell_max::Integer = 4, f_ref::Real, f_low::Real = 0.0)
    model in (:NRSur7dq4, :NRSur7dq4v2) || throw(ArgumentError("unsupported model $model"))
    2 <= ell_max <= (model === :NRSur7dq4 ? 4 : 5) || throw(ArgumentError("ell_max out of range"))
    f_ref >= f_low >= 0 || throw(ArgumentError("need f_ref >= f_low >= 0"))
    top = pyimport("gwsurrogate.surrogate")
    py_model = getproperty(top, String(model))(String(path))
    return PyCallGenerator(py_model, String(model), Int(ell_max), Float64(f_ref), Float64(f_low))
end

"""Return signed inertial modes in seconds and physical strain."""
function waveform_modes_td(g::PyCallGenerator, params::NamedTuple; dt = nothing, times = nothing)
    (dt === nothing || times === nothing) || throw(ArgumentError("specify at most one of dt and times"))
    M = params.total_mass
    tscale = MSUN_IN_SECONDS * M
    ascale = tscale * SPEED_OF_LIGHT / (1.0e6 * params.luminosity_distance * PARSEC)
    opts = PyDict()
    opts["init_orbphase"] = get(params, :init_orbphase, 0.0)
    haskey(params, :init_quat) && (opts["init_quat"] = collect(params.init_quat))
    kw = Dict{Symbol, Any}(
        :f_low => g.f_low * tscale, :f_ref => g.f_ref * tscale,
        :ellMax => g.ell_max, :precessing_opts => opts,
    )
    dt === nothing || (kw[:dt] = dt / tscale)
    times === nothing || (kw[:times] = collect(times) ./ tscale)
    q = params.mass_ratio
    spins = [collect(params.spin_1), collect(params.spin_2)]
    result = pycall(g.model, PyAny, q, spins[1], spins[2]; kw...)
    t = Vector{Float64}(result[1]) .* tscale
    hdict = result[2]
    labels = Tuple((l, m) for l in 2:g.ell_max for m in -l:l)
    modes = Matrix{ComplexF64}(undef, length(t), length(labels))
    for (j, label) in enumerate(labels)
        modes[:, j] .= Vector{ComplexF64}(hdict[label]) .* ascale
    end
    return (; t, mode_indices = labels, modes)
end

AbstractWaveforms.approximant(g::PyCallGenerator) = g.name
AbstractWaveforms.reference_frequency(g::PyCallGenerator) = g.f_ref
AbstractWaveforms.minimum_frequency(g::PyCallGenerator) = g.f_low
AbstractWaveforms.maximum_frequency(::PyCallGenerator) = nothing
AbstractWaveforms.spin_support(::PyCallGenerator) = :precessing
AbstractWaveforms.is_differentiable(::PyCallGenerator) = false
mode_angles(params, ::PyCallGenerator) = (params.inclination, oftype(float(params.phi_ref), π) / 2 - params.phi_ref)

end
