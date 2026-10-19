# SPDX-License-Identifier: GPL-3.0-or-later
# Derived from pySEOBNR 0.3.7; see ../NOTICE.md and ../provenance.toml.
# Modified: generic Julia arithmetic and one-based arrays.


function calc_isco_radius(a)
    # The exact Schwarzschild limit removes a removable 0/0 in its derivative.
    iszero(primal(a)) && return 6 - 4sqrt(2 / 3) * a
    z1 = 1 + (1 - a^2)^(1 / 3) * ((1 + a)^(1 / 3) + (1 - a)^(1 / 3))
    z2 = sqrt(3a^2 + z1^2)
    return 3 + z2 - sign(a) * sqrt((3 - z1) * (3 + z1 + 2z2))
end

function _bbh_HBR2016_ell(m1, m2, chi1z, chi2z, version)
    if (version == "M1J2")
        nJ = 2
        k00 = -(3.82116)
        k01 = -(1.2019)
        k02 = -(1.20764)
        k10 = 3.79245
        k11 = 1.18385
        k12 = 4.90494
        xi = 0.41616
    else
        if (version == "M3J3")
            nJ = 3
            k00 = -(5.83164)
            k01 = 2.87025
            k02 = -(1.53315)
            k03 = -(3.78893)
            k10 = 32.9127
            k11 = -(62.9901)
            k12 = 10.0068
            k13 = 56.1926
            k20 = -(136.832)
            k21 = 329.32
            k22 = -(13.2034)
            k23 = -(252.27)
            k30 = 210.075
            k31 = -(545.35)
            k32 = -(3.97509)
            k33 = 368.405
            xi = 0.463926
        else
            if (version == "M3J4")
                nJ = 4
                k00 = -(5.97723)
                k01 = 3.39221
                k02 = 4.48865
                k03 = -(5.77101)
                k04 = -(13.0459)
                k10 = 35.1278
                k11 = -(72.9336)
                k12 = -(86.0036)
                k13 = 93.7371
                k14 = 200.975
                k20 = -(146.822)
                k21 = 387.184
                k22 = 447.009
                k23 = -(467.383)
                k24 = -(884.339)
                k30 = 223.911
                k31 = -(648.502)
                k32 = -(697.177)
                k33 = 753.738
                k34 = 1166.89
                xi = 0.474046
            else
                throw(DomainError(r, "Hamiltonian domain"))
            end
        end
    end
    m = (m1 + m2)
    q = (m2 / m1)
    eta = ((m1 * m2) / (m * m))
    atot = ((chi1z + ((q * q) * chi2z)) / ((1.0 + q) * (1.0 + q)))
    aeff = (atot + ((xi * eta) * (chi1z + chi2z)))
    r_isco = calc_isco_radius(aeff)
    e_isco = ((1.0 - ((2.0 / 3.0) / r_isco))^0.5)
    l_isco = ((2.0 / (3.0 * (3.0^0.5))) * (1.0 + (2.0 * (((3.0 * r_isco) - 2.0)^0.5))))
    if (nJ >= 2)
        aeff2 = (aeff * aeff)
        ksum = (((k00 + (k01 * aeff)) + (k02 * aeff2)) + (eta * ((k10 + (k11 * aeff)) + (k12 * aeff2))))
    end
    if (nJ >= 3)
        eta2 = (eta * eta)
        eta3 = (eta2 * eta)
        aeff3 = (aeff2 * aeff)
        ksum = (((ksum + ((k03 + (eta * k13)) * aeff3)) + (eta2 * (((k20 + (k21 * aeff)) + (k22 * aeff2)) + (k23 * aeff3)))) + (eta3 * (((k30 + (k31 * aeff)) + (k32 * aeff2)) + (k33 * aeff3))))
    end
    if (nJ >= 4)
        ksum = (ksum + (((((k04 + (eta * k14)) + (eta2 * k24)) + (eta3 * k34)) * aeff3) * aeff))
    end
    ell = abs(((l_isco - ((2.0 * atot) * (e_isco - 1.0))) + (eta * ksum)))
    return ell
end

function bbh_final_spin_non_precessing_HBR2016(m1, m2, chi1z, chi2z, version = "M3J3")
    q = (m2 / m1)
    atot = ((chi1z + ((chi2z * q) * q)) / ((1.0 + q) * (1.0 + q)))
    ell = _bbh_HBR2016_ell(m1, m2, chi1z, chi2z, version)
    return (atot + (ell / (((1.0 / q) + 2.0) + q)))
end

function bbh_UIBfits_setup(m1, m2, chi1, chi2)
    m = (m1 + m2)
    msq = (m * m)
    m1sq = (m1 * m1)
    m2sq = (m2 * m2)
    eta = ((m1 * m2) / msq)
    eta2 = (eta * eta)
    eta3 = (eta2 * eta)
    eta4 = (eta2 * eta2)
    S1 = ((chi1 * m1sq) / msq)
    S2 = ((chi2 * m2sq) / msq)
    Stot = (S1 + S2)
    Shat = (((chi1 * m1sq) + (chi2 * m2sq)) / (m1sq + m2sq))
    Shat2 = (Shat * Shat)
    Shat3 = (Shat2 * Shat)
    Shat4 = (Shat2 * Shat2)
    chidiff = (chi1 - chi2)
    chidiff2 = (chidiff * chidiff)
    sqrt2 = (2.0^0.5)
    sqrt3 = (3.0^0.5)
    sqrt1m4eta = ((m1 - m2) / (m1 + m2))
    return (m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta)
end

function bbh_final_mass_non_precessing_UIB2016(m1, m2, chi1, chi2, version = "v2")
    (m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta) = bbh_UIBfits_setup(m1, m2, chi1, chi2)
    if (version == "v1")
        b10 = 0.487
        b20 = 0.295
        b30 = 0.17
        b50 = -(0.0717)
        a2 = 0.5635376058169299
        a3 = -(0.8661680065959881)
        a4 = 3.181941595301782
        b1 = -(0.15800074104558132)
        b2 = -(0.15815904609933157)
        b3 = -(0.14299315232521553)
        b5 = 8.908772171776285
        f20 = 3.8071100104582234
        f30 = 25.99956516423936
        f50 = 1.552929335555098
        f10 = 1.7004558922558886
        f21 = 0.0
        d10 = -(0.12282040108157262)
        d11 = -(3.499874245551208)
        d20 = 0.014200035799803777
        d30 = -(0.01873720734635449)
        d31 = -(5.1830734185518725)
        f11 = 14.39323998088354
        f31 = -(232.25752840151296)
        f51 = -(0.8427987782523847)
    else
        if (version == "v2")
            b10 = 0.346
            b20 = 0.211
            b30 = 0.128
            b50 = -(0.212)
            a2 = 0.5609904135313374
            a3 = -(0.84667563764404)
            a4 = 3.145145224278187
            b1 = -(0.2091189048177395)
            b2 = -(0.19709136361080587)
            b3 = -(0.1588185739358418)
            b5 = 2.9852925538232014
            f20 = 4.271313308472851
            f30 = 31.08987570280556
            f50 = 1.5673498395263061
            f10 = 1.8083565298668276
            f21 = 0.0
            d10 = -(0.09803730445895877)
            d11 = -(3.2283713377939134)
            d20 = 0.01118530335431078
            d30 = -(0.01978238971523653)
            d31 = -(4.91667749015812)
            f11 = 15.738082204419655
            f31 = -(243.6299258830685)
            f51 = -(0.5808669012986468)
        else
            throw(DomainError(r, "Hamiltonian domain"))
        end
    end
    Erad = ((((((((((1.0 + ((-(2.0) / 3.0) * sqrt2)) * eta) + (a2 * eta2)) + (a3 * eta3)) + (a4 * eta4)) * (((1.0 + (((b10 * b1) * Shat) * ((f10 + (f11 * eta)) + (((16.0 - (16.0 * f10)) - (4.0 * f11)) * eta2)))) + (((b20 * b2) * Shat2) * ((f20 + (f21 * eta)) + (((16.0 - (16.0 * f20)) - (4.0 * f21)) * eta2)))) + (((b30 * b3) * Shat3) * ((f30 + (f31 * eta)) + (((16.0 - (16.0 * f30)) - (4.0 * f31)) * eta2))))) / (1.0 + (((b50 * b5) * Shat) * ((f50 + (f51 * eta)) + (((16.0 - (16.0 * f50)) - (4.0 * f51)) * eta2))))) + ((((d10 * sqrt1m4eta) * eta2) * (1.0 + (d11 * eta))) * chidiff)) + (((((d30 * Shat) * sqrt1m4eta) * eta) * (1.0 + (d31 * eta))) * chidiff)) + ((d20 * eta3) * chidiff2))
    Mf = (m * (1.0 - Erad))
    return Mf
end
