"""
    hphc_fd(f::AbstractRange, params, generator::LALWaveformGenerator) -> (; plus, cross)

Frequency-domain plus and cross polarisations on the uniform grid `f` (Hz), following
bilby's `lal_binary_black_hole`. `params` is a `NamedTuple` with bilby's names and
units: required `mass_1`, `mass_2` (detector-frame solar masses),
`luminosity_distance` (Mpc), `theta_jn`, `phase`; optional (default 0) `a_1`, `a_2`,
`tilt_1`, `tilt_2`, `phi_12`, `phi_jl`, `eccentricity`.

`f` must start at 0 and have a positive step, as bilby assumes. The output is zero
outside `minimum_frequency ≤ f ≤ maximum_frequency`.

Throws [`LALError`](@ref) when LAL rejects the input, and `ArgumentError` for
malformed `params` or `f`; mass parameters other than `mass_1`, `mass_2` are converted
with `GWParameterConversion.convert_parameters`.
"""
function hphc_fd(f::AbstractRange, params, g::LALWaveformGenerator)
    p = bilby_parameters(params)
    check_frequency_grid(f)
    df = step(f)

    fns = lal()
    silent_errors!(fns)

    f_min, f_max = frequency_bounds(g, f)
    f_ref = g.reference_frequency
    approx = g.code

    start_frequency = if g.pn_amplitude_order != 0
        ccall(
            fns.XLALSimInspiralfLow2fStart, Float64, (Float64, Int32, Int32),
            f_min, g.pn_amplitude_order, approx,
        )
    else
        f_min
    end

    spins = bilby_to_lalsimulation_spins(p, f_ref)

    fn = g.fd_native ? fns.XLALSimInspiralChooseFDWaveform : fns.XLALSimInspiralFD
    hp = Ref{Ptr{COMPLEX16FrequencySeries}}(C_NULL)
    hc = Ref{Ptr{COMPLEX16FrequencySeries}}(C_NULL)
    try
        clear_errno!(fns)
        rc = ccall(
            fn, Cint,
            (
                Ptr{Ptr{COMPLEX16FrequencySeries}}, Ptr{Ptr{COMPLEX16FrequencySeries}},
                Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
                Float64, Float64, Float64, Float64, Float64, Float64,
                Float64, Float64, Float64, Float64,
                Ptr{Cvoid}, Cint,
            ),
            hp, hc,
            p.mass_1 * MSUN_SI, p.mass_2 * MSUN_SI,
            spins.s1x, spins.s1y, spins.s1z, spins.s2x, spins.s2y, spins.s2z,
            p.luminosity_distance * 1.0e6 * PC_SI, spins.iota, p.phase,
            0.0, p.eccentricity, 0.0,
            df, start_frequency, f_max, f_ref,
            g.dict.ptr, approx,
        )
        check(rc, g.fd_native ? :XLALSimInspiralChooseFDWaveform : :XLALSimInspiralFD, fns)

        plus = copy_masked(hp[], f, f_min, f_max)
        cross = copy_masked(hc[], f, f_min, f_max)
        if !g.fd_native
            # SimInspiralFD returns a series whose time origin is shifted; undo it as bilby does.
            ts = unsafe_load(hp[])
            dt = 1 / ts.deltaF + (ts.epoch.gpsSeconds + ts.epoch.gpsNanoSeconds * 1.0e-9)
            time_shift!(plus, f, f_min, f_max, dt)
            time_shift!(cross, f, f_min, f_max, dt)
        end
        return (; plus, cross)
    finally
        for h in (hp, hc)
            h[] == C_NULL || ccall(
                fns.XLALDestroyCOMPLEX16FrequencySeries, Cvoid,
                (Ptr{COMPLEX16FrequencySeries},), h[],
            )
        end
        # The generator owns the dictionary; keep it alive across the ccall.
        GC.@preserve g nothing
    end
end

# LAL may return fewer samples than the grid (zero-pad) or more, when padded to a
# power of two (truncate). Samples outside [f_min, f_max] are zeroed.
function copy_masked(series::Ptr{COMPLEX16FrequencySeries}, f, f_min, f_max)
    ts = unsafe_load(series)
    seq = unsafe_load(ts.data)
    out = zeros(ComplexF64, length(f))
    n = min(length(out), Int(seq.length))
    src = unsafe_wrap(Array, seq.data, Int(seq.length))
    @inbounds for i in 1:n
        fi = f[i]
        if f_min <= fi <= f_max
            out[i] = src[i]
        end
    end
    return out
end

function time_shift!(h, f, f_min, f_max, dt)
    # Same operation order as bilby: ((2π) * dt) * f.
    @inbounds for i in eachindex(h)
        fi = f[i]
        if f_min <= fi <= f_max
            h[i] *= cis(-(((2π) * dt) * fi))
        end
    end
    return h
end

"""
    validate_parameters(f::AbstractRange, params::NamedTuple, generator::LALWaveformGenerator) -> nothing

Check the grid and the parameter names and redundant mass parameters, as [`hphc_fd`](@ref) does
before calling LAL (which itself rejects unphysical values with a [`LALError`](@ref)). Provided
so generic code can validate with any backend.
"""
function validate_parameters(f::AbstractRange, params::NamedTuple, ::LALWaveformGenerator)
    check_frequency_grid(f)
    bilby_parameters(params)
    return nothing
end
