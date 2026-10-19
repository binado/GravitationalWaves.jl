"""Pinned pySEOBNR oracle. Run with `uv run --project python/seobnr python ...`.

No Julia runtime imports this environment. Large fixtures live in reference/seobnr.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import lal
import numpy as np
import tomlkit
from pyseobnr.generate_waveform import Ham_aligned_opt, generate_modes_opt
from pyseobnr.models.SEOBNRv5HM import SEOBNRv5HM_opt
from pyseobnr.eob.fits.fits_Hamiltonian import a6_NS, dSO, NR_deltaT, NR_deltaT_NS
from pyseobnr.eob.fits.GSF_fits import GSF_amplitude_fits
from pyseobnr.eob.fits.IV_fits import InputValueFits
from pyseobnr.eob.fits.MR_fits import MergerRingdownFits
from pyseobnr.eob.utils.nr_utils import bbh_final_mass_non_precessing_UIB2016, bbh_final_spin_non_precessing_HBR2016
from pyseobnr.eob.waveform.waveform import SEOBNRv5RRForce, compute_flux, compute_hlms
from pyseobnr.eob.dynamics.initial_conditions_aligned_opt import computeIC_opt

ROOT = Path(__file__).resolve().parents[2]
CASES = [(mass, q, s1, s2) for mass in (40, 60, 100)
         for q, s1, s2 in ((1, 0, 0), (2, .4, -.3), (4, -.6, .6))]


def verify_versions():
    for package, expected in (("pyseobnr", "0.3.7"), ("qnm", "0.4.4")):
        actual = importlib.metadata.version(package)
        if actual != expected:
            raise RuntimeError(f"Expected {package} {expected}, found {actual}")


def component_snapshot(q, s1, s2):
    omega0 = np.pi * 60 * lal.MTSUN_SI * 16
    model = SEOBNRv5HM_opt(q, s1, s2, omega0, Ham_aligned_opt, SEOBNRv5RRForce(),
                          settings={"postadiabatic": False, "return_modes": [(2, 2)]})
    model._initialize_params(phys_pars=model.phys_pars)
    model._set_H_coeffs()
    import re
    for key, value in GSF_amplitude_fits(model.nu).items():
        match = re.fullmatch(r"h(\d)(\d)_v(log)?(\d+)", key)
        l, m, log, power = match.groups()
        array = model.eob_pars.flux_params.extra_coeffs_log if log else model.eob_pars.flux_params.extra_coeffs
        array[int(l), int(m), int(power)] = value
    pos, mom = np.array([12., .3]), np.array([-.005, 4.])
    H = model.H
    args = (pos, mom, s1, s2, model.m_1, model.m_2)
    dynamics = H.dynamics(*args)
    omega_circ = H.omega(pos, np.array([0., 4.]), *args[2:])
    flux = compute_flux(*pos, *mom, dynamics[3], omega_circ, model.nu * dynamics[4], model.eob_pars)
    iv = InputValueFits(model.m_1, model.m_2, [0, 0, s1], [0, 0, s2])
    mr = MergerRingdownFits(model.m_1, model.m_2, [0, 0, s1], [0, 0, s2])
    result = dict(q=float(q), chi1=float(s1), chi2=float(s2),
                  a6=float(a6_NS(model.nu)), dSO=float(dSO(model.nu, model.ap, model.am)),
                  delta_t=float(NR_deltaT_NS(model.nu) + NR_deltaT(model.nu, model.ap, model.am)),
                  hamiltonian=list(H(*args, verbose=True)), gradient=list(H.grad(*args)),
                  hessian=np.asarray(H.hessian(*args)).reshape(-1).tolist(),
                  dynamics=list(dynamics), omega_circ=float(omega_circ), flux=float(flux),
                  final_mass=float(bbh_final_mass_non_precessing_UIB2016(model.m_1, model.m_2, s1, s2)),
                  final_spin=float(bbh_final_spin_non_precessing_HBR2016(model.m_1, model.m_2, s1, s2, version="M3J4")))
    for name in ("habs", "hdot", "hdotdot", "omega", "omegadot"):
        result["iv_" + name] = float(getattr(iv, name)()[2, 2])
    for name in ("c1f", "c2f", "d1f", "d2f"):
        result["mr_" + name] = float(getattr(mr, name)()[2, 2])
    mode = np.array(compute_hlms(np.array([[*pos, *mom, dynamics[4], dynamics[3], omega_circ]]), model.eob_pars)[2, 2])[0]
    result["mode_real"], result["mode_imag"] = float(mode.real), float(mode.imag)
    result["initial_conditions"] = list(computeIC_opt(omega0, H, model.RR, s1, s2, model.m_1, model.m_2, params=model.eob_pars))
    result["initial_omega"] = float(omega0)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--components", action="store_true")
    parser.add_argument("--waves", action="store_true")
    args = parser.parse_args()
    verify_versions()
    if args.components:
        data = {"pyseobnr_version": "0.3.7", "qnm_version": "0.4.4",
                "cases": [component_snapshot(q, s1, s2) for _, q, s1, s2 in CASES[:3]]}
        path = ROOT / "packages/GWSEOBNR/test/oracle.toml"
        path.write_text(tomlkit.dumps(data))
        print(path, flush=True)
    if args.waves:
        dest = ROOT / "reference/seobnr"
        dest.mkdir(parents=True, exist_ok=True)
        for mass, q, s1, s2 in CASES:
            settings = {"M": mass, "dt": 1/4096, "f_ref": 20,
                        "postadiabatic": False, "return_modes": [(2, 2)],
                        "rtol_ode": 1e-11, "atol_ode": 1e-12}
            t, modes, model = generate_modes_opt(q, s1, s2, np.pi*mass*lal.MTSUN_SI*16,
                                               settings=settings, debug=True)
            name = f"M{mass}_q{q}"
            np.savez(dest / (name + ".npz"), t=t, h22=modes["2,2"], dynamics=model.dynamics,
                     t_isco=model.t_ISCO, t_attach=model.t_attach, r_isco=model.r_ISCO,
                     attachment_check=model.attachment_check, nqc=np.array(list(model.nqc_coeffs[2, 2].values())))
            print(name, "attachment fallback:", model.attachment_check, flush=True)


if __name__ == "__main__":
    main()
