# SPDX-License-Identifier: GPL-3.0-or-later
# Ported from pySEOBNR 0.3.7 waveform.pyx and GSF_fits.py; see NOTICE.md.
# Modified to precompute draw-dependent coefficients in their promoted scalar type.
function multipole_prefix(m1, m2, l, m)
    M = m1 + m2
    x1, x2 = m1 / M, m2 / M
    eta = m1 * m2 / (M * M)
    epsilon = (l + m) % 2
    s = iseven(m) ? 1 : -1
    c = if m1 != m2 || s == 1
        x2^(l + epsilon - 1) + s * x1^(l + epsilon - 1)
    else
        l <= 3 ? -1.0 : l <= 5 ? -0.5 : 0.0
    end
    n = (im * m)^l
    dfac = prod(1:2:(2l + 1))
    if epsilon == 0
        n *= 8.0 * π / dfac
        n *= sqrt(((l + 1) * (l + 2)) / (l * (l - 1)))
    else
        n = -n
        n *= im * 16.0 * π / dfac
        n *= sqrt(((2l + 1) * (l + 2) * (l * l - m * m)) / ((2l - 1) * (l + 1) * l * (l - 1)))
    end
    return n * eta * c
end

function flux_coefficients(c)
    T = typeof(c.nu)
    rho, rholog, f, fvh = (zeros(T, 9, 9, 11) for _ in 1:4)
    delta, deltavh = (zeros(T, 9, 9, 11) for _ in 1:2)
    chiS, chiA = (c.chi_1 + c.chi_2) / 2, (c.chi_1 - c.chi_2) / 2
    dm = c.X_1 - c.X_2
    a = (1 - 2c.nu) * chiS + dm * chiA
    rho_coefficients!(c.nu, dm, a, chiS, chiA, rho, rholog, f, fvh, true)
    delta_coefficients!(c.nu, dm, a, chiS, chiA, delta, deltavh)
    # Nonspinning second-order self-force amplitude corrections.
    gsf = (
        (2, 2, 8, 21.2), (2, 2, 10, -411.0),
        (3, 3, 8, 12.0), (3, 3, 10, -215.0),
        (2, 1, 6, 1.65), (2, 1, 8, 26.5), (2, 1, 10, 80.0),
        (4, 4, 6, -3.56), (4, 4, 8, 15.6), (4, 4, 10, -216.0),
        (5, 5, 4, -2.61), (5, 5, 6, 1.25), (5, 5, 8, -35.7),
        (3, 2, 6, 0.333), (3, 2, 8, -6.5),
        (3, 2, 10, 98 - (1312549797426453052 / 176264081083715625) / c.nu),
        (4, 3, 4, -0.654), (4, 3, 6, -3.69),
        (4, 3, 8, 18.5 - (2465107182496333 / 460490801971200) / c.nu),
    )
    extra, extralog = zeros(T, 9, 9, 11), zeros(T, 9, 9, 11)
    for (l, m, v, coefficient) in gsf
        extra[l + 1, m + 1, v + 1] = c.nu * coefficient
    end
    extralog[4, 3, 11] = c.nu * ((18778864 / 12629925) / c.nu)
    extralog[5, 4, 9] = c.nu * ((174381 / 67760) / c.nu)
    prefixes = zeros(Complex{T}, 9, 9)
    for l in 2:8, m in 1:l
        prefixes[l + 1, m + 1] = multipole_prefix(c.X_1, c.X_2, l, m)
    end
    return (; rho, rholog, f, fvh, delta, deltavh, extra, extralog, prefixes)
end

function residual_amplitude(vs, vh, l, m, c, coeffs)
    eulerlog = EULER_GAMMA + log(2.0 * m * vs[2])
    rho, f = one(vh), zero(vh)
    for j in 1:10
        rho += (
            coeffs.rho[l + 1, m + 1, j + 1] + coeffs.extra[l + 1, m + 1, j + 1] +
                (coeffs.rholog[l + 1, m + 1, j + 1] + coeffs.extralog[l + 1, m + 1, j + 1]) * eulerlog
        ) * vs[j + 1]
        if isodd(m)
            f += coeffs.f[l + 1, m + 1, j + 1] * vs[j + 1]
        end
    end
    extra = l == 3 && m == 3 ? (vh^3)^2 * coeffs.fvh[4, 4, 7] * im : zero(complex(vh))
    ff = f + extra
    return abs(c.nu - 0.25) < 1.0e-14 && isodd(m) ? ff : rho^l + ff
end

"""Full upstream radiation flux, including every 2 ≤ ℓ ≤ 8 multipole."""
function radiation_flux(pphi, omega, omega_circ, H, c, coeffs)
    v = omega^(1 / 3)
    vh = (H * omega)^(1 / 3)
    vs = ntuple(j -> v^(j - 1), 11)
    source1 = (H * H - 1.0) / (2.0 * c.nu) + 1.0
    source2 = v * pphi
    v_phi = omega / omega_circ^(2 / 3)
    v_phi2 = v_phi * v_phi
    flux = zero(v)
    for m in 1:8
        khat = H * m * omega
        khat4 = 4.0 * khat * khat
        fourpi = 4.0 * π * khat
        prefac = sqrt(fourpi / (1.0 - exp(-fourpi)))
        tailprod = one(v)
        for l in 1:8
            tailprod *= khat4 + l^2
            l >= max(m, 2) || continue
            tail = prefac * sqrt(tailprod) / factorial(l)
            epsilon = (l + m) % 2
            newton2 = abs2(coeffs.prefixes[l + 1, m + 1]) * v_phi2^(l + epsilon) * YLMS[m + 1, l - epsilon + 1]^2
            source = epsilon == 0 ? source1 : source2
            residual = residual_amplitude(vs, vh, l, m, c, coeffs)
            flux += m * m * omega * omega * tail^2 * source^2 * newton2 * abs2(residual)
        end
    end
    return flux / (8π)
end

function radiation_reaction(pr, pphi, omega, omega_circ, H, c, coeffs)
    flux = radiation_flux(pphi, omega, omega_circ, c.nu * H, c, coeffs) / c.nu
    f_over_om = flux / omega
    return -pr / pphi * f_over_om, -f_over_om
end
