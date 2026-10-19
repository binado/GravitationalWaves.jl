function qnm_polynomial(a, knots, coeffs)
    0 <= primal(a) <= knots[end] || throw(ArgumentError("remnant spin outside bundled QNM interpolation domain"))
    i = clamp(searchsortedlast(knots, primal(a)), 1, length(knots) - 1)
    x = a - knots[i]
    return ((coeffs[i, 1] * x + coeffs[i, 2]) * x + coeffs[i, 3]) * x + coeffs[i, 4]
end

"""Pinned qnm 0.4.4 (2,2,0) interpolant, for either sign of remnant spin."""
function qnm_220(af, Mf)
    if af > 0
        re = qnm_polynomial(af, QNM_POSITIVE_REAL_KNOTS, QNM_POSITIVE_REAL_COEFFS)
        im = qnm_polynomial(af, QNM_POSITIVE_IMAG_KNOTS, QNM_POSITIVE_IMAG_COEFFS)
    else
        re = qnm_polynomial(-af, QNM_NEGATIVE_REAL_KNOTS, QNM_NEGATIVE_REAL_COEFFS)
        im = qnm_polynomial(-af, QNM_NEGATIVE_IMAG_KNOTS, QNM_NEGATIVE_IMAG_COEFFS)
    end
    return complex(re, im) / Mf
end
