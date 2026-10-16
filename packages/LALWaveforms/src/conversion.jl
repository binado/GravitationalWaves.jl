"""
    bilby_to_lalsimulation_spins(p::BBHParameters, reference_frequency)

Port of bilby's `bilby_to_lalsimulation_spins` (`bilby/gw/conversion.py`). Returns
`(; iota, s1x, s1y, s1z, s2x, s2y, s2z)`; the spins are dimensionless Cartesian
components in the LAL frame at `reference_frequency`. The aligned-spin short circuit is
shared with GWParameterConversion; the precessing branch calls LAL itself, which keeps it the
oracle for `GWParameterConversion.transform_precessing_initial_conditions`.
"""
function bilby_to_lalsimulation_spins(p::BBHParameters, reference_frequency::Real)
    if is_aligned(p.a_1, p.tilt_1) && is_aligned(p.a_2, p.tilt_2)
        return aligned_spins_to_lal_frame(p.theta_jn, p.a_1, p.tilt_1, p.a_2, p.tilt_2)
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
