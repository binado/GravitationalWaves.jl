# `XLALSimInspiralfLow2fStart`: the amplitude corrections of
# TaylorF2 start earlier than the requested frequency.
function f_low_to_f_start(f_low, amplitude_order::Integer)
    ampOrder = amplitude_order == -1 ? 6 : amplitude_order # MAX_NONPRECESSING_AMP_PN_ORDER
    return f_low * 2.0 / (ampOrder + 2)
end

# All parameters converted to one floating-point type, so the kernels are type stable and dual
# numbers in any parameter propagate through all of them.
function promote_parameters(p::NamedTuple)
    T = float(promote_type(map(typeof, values(p))...))
    return T, map(T, p)
end

"""
    hphc_fd(f::AbstractRange, params::NamedTuple, settings::WaveformSettings) -> (; plus, cross)

Frequency-domain plus and cross polarisations on the uniform grid `f` (Hz), following bilby's
`lal_binary_black_hole` and the interface of `LALWaveforms.hphc_fd`. `params` has bilby's
names and units: required `mass_1`, `mass_2` (detector-frame solar masses; other mass
parameters are converted by `GWParameterConversion.convert_parameters`),
`luminosity_distance` (Mpc), `theta_jn`, `phase`; optional (default 0) `a_1`, `a_2`, `tilt_1`,
`tilt_2`, `phi_12`, `phi_jl`, `eccentricity`.

`f` must start at 0 and have a positive step, as bilby assumes. The output is zero outside
`minimum_frequency ≤ f ≤ maximum_frequency`. The spins must be (anti-)aligned with the orbital
angular momentum and `eccentricity` zero, else `ArgumentError` is thrown.

The element type is `Complex{T}` with `T` the floating-point type promoted from the
parameters, so dual numbers (ForwardDiff) and Enzyme or Mooncake reverse mode pass through.
"""
function hphc_fd(f::AbstractRange, params::NamedTuple, settings::WaveformSettings)
    p0 = aligned_parameters(params)
    first(f) == 0 || throw(ArgumentError("frequency grid must start at 0, got $(first(f))"))
    df = step(f)
    df > 0 || throw(ArgumentError("frequency grid must have a positive step, got $df"))
    length(f) >= 2 || throw(ArgumentError("frequency grid needs at least two points"))
    T, p = promote_parameters(p0)

    f_min = settings.minimum_frequency
    f_max = something(settings.maximum_frequency, Float64(last(f)))
    f_start = settings.pn_amplitude_order != 0 ?
        f_low_to_f_start(f_min, settings.pn_amplitude_order) : f_min

    h = waveform_fd(settings.approximant, T, length(f), Float64(df), p, settings, f_start, f_max)

    # The models return h(f) for ι = 0; other orientations multiply by the inclination
    # dependence, hc ∝ -i hp (ChooseFDWaveform in LALSimInspiral).
    cfac = cos(p.inclination)
    pfac = 0.5 * (1.0 + cfac * cfac)
    plus = similar(h)
    cross = similar(h)
    for i in eachindex(h)
        fi = f[i]
        if f_min <= fi <= f_max
            x = h[i]
            plus[i] = complex(pfac * real(x), pfac * imag(x))
            cross[i] = complex(cfac * imag(x), -(cfac * real(x)))
        else
            plus[i] = zero(eltype(h))
            cross[i] = zero(eltype(h))
        end
    end
    return (; plus, cross)
end

waveform_fd(::TaylorF2, args...) = taylorf2_fd(args...)
