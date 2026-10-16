"""
    bilby_to_lalsimulation_spins(p::BBHParameters, reference_frequency)

Port of bilby's `bilby_to_lalsimulation_spins` (`bilby/gw/conversion.py`). Returns
`(; iota, s1x, s1y, s1z, s2x, s2y, s2z)`; the spins are dimensionless Cartesian
components in the LAL frame at `reference_frequency`.
"""
function bilby_to_lalsimulation_spins(p::BBHParameters, reference_frequency::Real)
    # Compare with the rounded double (np.pi), not the Irrational: `Float64(π) == π` is false.
    pi64 = Float64(π)
    aligned1 = p.a_1 == 0 || p.tilt_1 == 0 || p.tilt_1 == pi64
    aligned2 = p.a_2 == 0 || p.tilt_2 == 0 || p.tilt_2 == pi64
    if aligned1 && aligned2
        return (;
            iota = p.theta_jn, s1x = 0.0, s1y = 0.0, s1z = p.a_1 * cos(p.tilt_1),
            s2x = 0.0, s2y = 0.0, s2z = p.a_2 * cos(p.tilt_2),
        )
    end
    return transform_precessing(p, Float64(reference_frequency))
end

function transform_precessing(p::BBHParameters, reference_frequency::Float64)
    fns = lal()
    silent_errors!(fns)
    iota, s1x, s1y, s1z, s2x, s2y, s2z = ntuple(_ -> Ref(0.0), 7)
    clear_errno!(fns)
    rc = ccall(
        fns.XLALSimInspiralTransformPrecessingNewInitialConditions, Cint,
        (
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64},
            Float64, Float64, Float64, Float64, Float64, Float64, Float64,
            Float64, Float64, Float64, Float64,
        ),
        iota, s1x, s1y, s1z, s2x, s2y, s2z,
        p.theta_jn, p.phi_jl, p.tilt_1, p.tilt_2, p.phi_12, p.a_1, p.a_2,
        p.mass_1 * MSUN_SI, p.mass_2 * MSUN_SI, reference_frequency, p.phase,
    )
    check(rc, :XLALSimInspiralTransformPrecessingNewInitialConditions, fns)
    return (;
        iota = iota[], s1x = s1x[], s1y = s1y[], s1z = s1z[],
        s2x = s2x[], s2y = s2y[], s2z = s2z[],
    )
end
