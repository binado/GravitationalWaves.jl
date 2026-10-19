# SPDX-License-Identifier: GPL-3.0-or-later
# Derived from pySEOBNR 0.3.7; see ../NOTICE.md and ../provenance.toml.
# Modified: generic Julia arithmetic and one-based arrays.

function dSO(nu, ap, am)
    return ((((((((((((((((((((-(7.71251231383957) * (am^3)) - ((17.2294679794015 * (am^2)) * ap)) - ((238.430383378296 * (am^2)) * nu)) + (69.5461667822545 * (am^2))) - ((10.5225438990315 * am) * (ap^2))) + (((362.767393298729 * am) * ap) * nu)) - ((85.8036338010274 * am) * ap)) - ((1254.66845939312 * am) * (nu^2))) + ((472.431937787377 * am) * nu)) - (39.742317057316 * am)) - (7.58458103577458 * (ap^3))) - ((42.7601129678844 * (ap^2)) * nu)) + (18.1783435552183 * (ap^2))) - ((201.905934468847 * ap) * (nu^2))) - ((90.5790079104259 * ap) * nu)) + (49.6299175121658 * ap)) + (478.546231305475 * (nu^3))) + (679.521769948995 * (nu^2))) - (177.334831768076 * nu)) - 37.6897780220529)
end

function NR_deltaT(nu, ap, am)
    return ((nu^((-(1.0) / 5) + (0 * nu))) * (((((((((((((((((8.39238879807543 * (am^2)) * ap) - ((16.9056858928167 * (am^2)) * nu)) + (7.23410583477034 * (am^2))) + ((6.38975598319936 * am) * (ap^2))) + (((179.569824846781 * am) * ap) * nu)) - ((40.6063653476775 * am) * ap)) + ((144.253395844761 * am) * (nu^2))) - ((90.1929138487509 * am) * nu)) + (14.2203101910927 * am)) - (6.78913884987037 * (ap^4))) + (5.39962303470497 * (ap^3))) - ((132.224950777226 * (ap^2)) * nu)) + (49.8016443361381 * (ap^2))) + ((384.201018794943 * ap) * (nu^2))) - ((141.253181790353 * ap) * nu)) + (17.5710132409988 * ap)))
end

function a6_NS(nu)
    par = [41.7877875, -(3021.93382), 33414.4394, -(169019.14), 329523.262]
    return ((((par[0 + 1] + (par[1 + 1] * nu)) + (par[2 + 1] * (nu^2))) + (par[3 + 1] * (nu^3))) + (par[4 + 1] * (nu^4)))
end

function NR_deltaT_NS(nu)
    par = [10.0513217, -(59.62318), -(1056.87385), -(9793.17619), 55565.2392]
    return ((nu^((-(1.0) / 5) + (par[0 + 1] * nu))) * (((par[1 + 1] + (par[2 + 1] * nu)) + (par[3 + 1] * (nu^2))) + (par[4 + 1] * (nu^3))))
end

"""Draw-dependent aligned Hamiltonian coefficients; masses are normalized to total mass."""
function hamiltonian_coefficients(m1, m2, chi1, chi2)
    X_1, X_2, chi_1, chi_2 = promote(m1 / (m1 + m2), m2 / (m1 + m2), chi1, chi2)
    nu = X_1 * X_2
    ap, am = X_1 * chi_1 + X_2 * chi_2, X_1 * chi_1 - X_2 * chi_2
    return (; X_1, X_2, chi_1, chi_2, nu, a6 = a6_NS(nu), dSO = dSO(nu, ap, am))
end
