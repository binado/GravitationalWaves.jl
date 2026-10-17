# Owning wrapper for a `LALDict*`. The finalizer frees it once no generator
# refers to it any more.
mutable struct LALDictHandle
    ptr::Ptr{Cvoid}
    function LALDictHandle(ptr::Ptr{Cvoid})
        h = new(ptr)
        finalizer(h) do h
            if h.ptr != C_NULL
                ccall(lal().XLALDestroyDict, Cvoid, (Ptr{Cvoid},), h.ptr)
                h.ptr = C_NULL
            end
        end
        return h
    end
end

"""
    LALWaveformGenerator(approximant; kwargs...)
    LALWaveformGenerator(; approximant, kwargs...)

Everything about a waveform call that does not vary between parameter draws: the LAL
approximant, frequency bounds and the LAL waveform dictionary (post-Newtonian orders, mode
array). Build it once and reuse it; it is immutable and its dictionary is only read after
construction. `approximant` is any name LAL knows (a `Symbol` also works), else [`LALError`](@ref)
is thrown.

Defaults are those of bilby's `lal_binary_black_hole`.

# Keywords
- `reference_frequency = 50.0`, `minimum_frequency = 20.0` (Hz)
- `maximum_frequency = nothing`: the last frequency of the grid
- `pn_spin_order = -1`, `pn_tidal_order = -1`, `pn_phase_order = -1`,
  `pn_amplitude_order = 0`
- `mode_array = nothing`: vector of `(l, m)` tuples restricting the modes used
"""
struct LALWaveformGenerator <: AbstractWaveformGenerator
    name::String
    code::Cint
    fd_native::Bool
    spin_support::Symbol
    reference_frequency::Float64
    minimum_frequency::Float64
    maximum_frequency::Union{Nothing, Float64}
    pn_amplitude_order::Int
    dict::LALDictHandle
end

function LALWaveformGenerator(
        approximant::Union{Symbol, AbstractString};
        reference_frequency::Real = 50.0,
        minimum_frequency::Real = 20.0,
        maximum_frequency::Union{Nothing, Real} = nothing,
        pn_spin_order::Integer = -1,
        pn_tidal_order::Integer = -1,
        pn_phase_order::Integer = -1,
        pn_amplitude_order::Integer = 0,
        mode_array = nothing,
    )
    fns = lal()
    silent_errors!(fns)
    name = String(approximant)

    clear_errno!(fns)
    code = ccall(
        fns.XLALSimInspiralGetApproximantFromString, Cint, (Cstring,), name
    )
    check(code < 0 ? -1 : 0, :XLALSimInspiralGetApproximantFromString, fns)
    fd_native = ccall(
        fns.XLALSimInspiralImplementedFDApproximants, Cint, (Cint,), code
    ) != 0
    spin = ccall(fns.XLALSimInspiralGetSpinSupportFromApproximant, Cint, (Cint,), code)

    dict = new_dict(
        fns; pn_spin_order, pn_tidal_order, pn_phase_order, pn_amplitude_order, mode_array
    )
    return LALWaveformGenerator(
        name, code, fd_native, spin_support_symbol(spin),
        reference_frequency, minimum_frequency,
        maximum_frequency === nothing ? nothing : Float64(maximum_frequency),
        pn_amplitude_order, dict,
    )
end

LALWaveformGenerator(; approximant, kwargs...) = LALWaveformGenerator(approximant; kwargs...)

# `SimInspiralSpinSupport` in LALSimInspiral.h: SPINLESS, SINGLESPIN, ALIGNEDSPIN,
# PRECESSINGSPIN, CASEBYCASE_SPINSUPPORT. Single-spin and case-by-case models are not
# restricted to aligned spins, so they count as precessing.
function spin_support_symbol(code::Integer)
    return code == 0 ? :none : code == 2 ? :aligned : :precessing
end

approximant(g::LALWaveformGenerator) = g.name
reference_frequency(g::LALWaveformGenerator) = g.reference_frequency
minimum_frequency(g::LALWaveformGenerator) = g.minimum_frequency
maximum_frequency(g::LALWaveformGenerator) = g.maximum_frequency
spin_support(g::LALWaveformGenerator) = g.spin_support
supports_mode_array(::LALWaveformGenerator) = true
is_differentiable(::LALWaveformGenerator) = false

# Mirrors bilby's `set_waveform_dictionary` for the options we support: the four PN
# orders, zero tidal deformabilities, and an optional mode array.
function new_dict(
        fns::LALFunctions; pn_spin_order, pn_tidal_order, pn_phase_order,
        pn_amplitude_order, mode_array,
    )
    ptr = ccall(fns.XLALCreateDict, Ptr{Cvoid}, ())
    ptr == C_NULL && throw(LALError(:XLALCreateDict, 0, "could not allocate LALDict"))
    handle = LALDictHandle(ptr)

    for (sym, value) in (
            (fns.XLALSimInspiralWaveformParamsInsertPNSpinOrder, pn_spin_order),
            (fns.XLALSimInspiralWaveformParamsInsertPNTidalOrder, pn_tidal_order),
            (fns.XLALSimInspiralWaveformParamsInsertPNPhaseOrder, pn_phase_order),
            (fns.XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder, pn_amplitude_order),
        )
        clear_errno!(fns)
        rc = ccall(sym, Cint, (Ptr{Cvoid}, Int32), ptr, value)
        check(rc, :XLALSimInspiralWaveformParamsInsertPNOrder, fns)
    end
    for sym in (
            fns.XLALSimInspiralWaveformParamsInsertTidalLambda1,
            fns.XLALSimInspiralWaveformParamsInsertTidalLambda2,
        )
        clear_errno!(fns)
        rc = ccall(sym, Cint, (Ptr{Cvoid}, Float64), ptr, 0.0)
        check(rc, :XLALSimInspiralWaveformParamsInsertTidalLambda, fns)
    end

    if mode_array !== nothing
        modes = ccall(fns.XLALSimInspiralCreateModeArray, Ptr{Cvoid}, ())
        modes == C_NULL && throw(LALError(:XLALSimInspiralCreateModeArray, 0, "allocation failed"))
        try
            for (l, m) in mode_array
                clear_errno!(fns)
                r = ccall(
                    fns.XLALSimInspiralModeArrayActivateMode, Ptr{Cvoid},
                    (Ptr{Cvoid}, Cuint, Cint), modes, l, m,
                )
                check(r == C_NULL ? -1 : 0, :XLALSimInspiralModeArrayActivateMode, fns)
            end
            clear_errno!(fns)
            rc = ccall(
                fns.XLALSimInspiralWaveformParamsInsertModeArray, Cint,
                (Ptr{Cvoid}, Ptr{Cvoid}), ptr, modes,
            )
            check(rc, :XLALSimInspiralWaveformParamsInsertModeArray, fns)
        finally
            # The dictionary stores its own copy of the value.
            ccall(fns.XLALDestroyValue, Cvoid, (Ptr{Cvoid},), modes)
        end
    end
    return handle
end
