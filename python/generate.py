"""Write bilby reference waveforms for every case in cases.json to test/reference/*.npz."""

from __future__ import annotations

import json
from pathlib import Path

import bilby.gw.conversion as conv
import bilby.gw.source as source
import numpy as np
from bilby.core.utils import solar_mass

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "reference"

# (mass_1, mass_2) detector-frame pairs for the mass-conversion fixtures; m1 >= m2.
MASS_PAIRS = [(36.0, 29.0), (50.0, 20.0), (30.0, 30.0), (10.5, 1.4), (80.3, 12.1), (1.4, 1.35)]


# Spin-frame cases: bilby parameter names, reference frequency in Hz.
SPIN_CASES = [
    {"mass_1": 36.0, "mass_2": 29.0, "theta_jn": 0.4, "phase": 1.3, "a_1": 0.3, "a_2": 0.2, "tilt_1": 0.0, "tilt_2": float(np.pi), "f_ref": 50.0},
    {"mass_1": 40.0, "mass_2": 25.0, "theta_jn": 0.9, "phase": 2.1, "a_1": 0.6, "a_2": 0.4, "tilt_1": 0.7, "tilt_2": 1.9, "phi_12": 1.1, "phi_jl": 2.5, "f_ref": 30.0},
    {"mass_1": 50.0, "mass_2": 20.0, "theta_jn": 1.2, "phase": 0.5, "a_1": 0.5, "a_2": 0.3, "tilt_1": 1.0, "tilt_2": 0.6, "phi_12": 0.3, "phi_jl": 4.0, "f_ref": 50.0},
    {"mass_1": 10.5, "mass_2": 1.4, "theta_jn": 2.6, "phase": 5.0, "a_1": 0.9, "a_2": 0.05, "tilt_1": 2.9, "tilt_2": 0.2, "phi_12": 5.8, "phi_jl": 0.3, "f_ref": 20.0},
    {"mass_1": 30.0, "mass_2": 30.0, "theta_jn": 0.05, "phase": 3.0, "a_1": 0.95, "a_2": 0.95, "tilt_1": 1.5, "tilt_2": 1.6, "phi_12": 3.1, "phi_jl": 6.0, "f_ref": 100.0},
    {"mass_1": 80.3, "mass_2": 12.1, "theta_jn": 3.0, "phase": 0.0, "a_1": 0.2, "a_2": 0.7, "tilt_1": 0.4, "tilt_2": 2.4, "phi_12": 2.0, "phi_jl": 1.0, "f_ref": 40.0},
]


def write_spin_cases() -> list[dict]:
    """bilby's `bilby_to_lalsimulation_spins` for each case in SPIN_CASES."""
    out = []
    for case in SPIN_CASES:
        c = {"phi_12": 0.0, "phi_jl": 0.0} | case
        spins = conv.bilby_to_lalsimulation_spins(
            theta_jn=c["theta_jn"],
            phi_jl=c["phi_jl"],
            tilt_1=c["tilt_1"],
            tilt_2=c["tilt_2"],
            phi_12=c["phi_12"],
            a_1=c["a_1"],
            a_2=c["a_2"],
            mass_1=c["mass_1"] * solar_mass,
            mass_2=c["mass_2"] * solar_mass,
            reference_frequency=c["f_ref"],
            phase=c["phase"],
        )
        out.append({**c, "spins": [float(x) for x in spins]})  # iota, s1x, s1y, s1z, s2x, s2y, s2z
    return out


def finite_or_none(x) -> float | None:
    """bilby returns nan (numpy) or a complex number (Python floats) where its closed-form cubic
    root has no real solution; the fixtures record those cases as null."""
    x = complex(x)
    return float(x.real) if x.imag == 0 and np.isfinite(x.real) else None


def write_conversion_cases() -> None:
    """Write bilby mass-conversion results to reference/conversion.json.

    Inverse functions are fed bilby's own forward outputs so each one is tested in isolation.
    """
    cases = []
    for m1, m2 in MASS_PAIRS:
        mc = float(conv.component_masses_to_chirp_mass(m1, m2))
        mt = float(conv.component_masses_to_total_mass(m1, m2))
        q = float(conv.component_masses_to_mass_ratio(m1, m2))
        eta = float(conv.component_masses_to_symmetric_mass_ratio(m1, m2))
        cases.append(
            {
                "mass_1": m1,
                "mass_2": m2,
                "chirp_mass": mc,
                "total_mass": mt,
                "mass_ratio": q,
                "symmetric_mass_ratio": eta,
                "from_total_mass_mass_ratio": [
                    float(x) for x in conv.total_mass_and_mass_ratio_to_component_masses(q, mt)
                ],
                "from_chirp_mass_mass_ratio": [
                    float(x) for x in conv.chirp_mass_and_mass_ratio_to_component_masses(mc, q)
                ],
                "total_mass_from_chirp_mass_mass_ratio": float(
                    conv.chirp_mass_and_mass_ratio_to_total_mass(mc, q)
                ),
                "eta_from_chirp_mass_total_mass": float(
                    conv.chirp_mass_and_total_mass_to_symmetric_mass_ratio(mc, mt)
                ),
                "mass_ratio_from_eta": float(conv.symmetric_mass_ratio_to_mass_ratio(eta)),
                "mass_ratio_from_mass_1_chirp_mass": finite_or_none(
                    conv.mass_1_and_chirp_mass_to_mass_ratio(np.float64(m1), np.float64(mc))
                ),
                "mass_ratio_from_mass_2_chirp_mass": finite_or_none(
                    conv.mass_2_and_chirp_mass_to_mass_ratio(np.float64(m2), np.float64(mc))
                ),
            }
        )
    OUT.mkdir(exist_ok=True)
    (OUT / "conversion.json").write_text(json.dumps({"masses": cases, "spins": write_spin_cases()}, indent=1))
    print(f"wrote conversion.json: {len(cases)} mass cases, {len(SPIN_CASES)} spin cases")


def main() -> None:
    spec = json.loads((HERE / "cases.json").read_text())
    df = spec["frequency"]["delta_f"]
    f_max = spec["frequency"]["f_max"]
    frequency = np.arange(int(round(f_max / df)) + 1) * df

    OUT.mkdir(exist_ok=True)
    write_conversion_cases()
    for case in spec["cases"]:
        # bilby's signature has no defaults for the spin parameters; the Julia side defaults them to 0.
        zeros = dict.fromkeys(["a_1", "a_2", "tilt_1", "tilt_2", "phi_12", "phi_jl"], 0.0)
        params, settings = zeros | case["params"], dict(case["settings"])
        if "mode_array" in settings:
            settings["mode_array"] = [tuple(m) for m in settings["mode_array"]]
        h = source.lal_binary_black_hole(
            frequency,
            **params,
            waveform_approximant=case["approximant"],
            **settings,
        )
        spins = conv.bilby_to_lalsimulation_spins(
            theta_jn=params["theta_jn"],
            phi_jl=params.get("phi_jl", 0.0),
            tilt_1=params.get("tilt_1", 0.0),
            tilt_2=params.get("tilt_2", 0.0),
            phi_12=params.get("phi_12", 0.0),
            a_1=params.get("a_1", 0.0),
            a_2=params.get("a_2", 0.0),
            mass_1=params["mass_1"] * solar_mass,
            mass_2=params["mass_2"] * solar_mass,
            reference_frequency=settings.get("reference_frequency", 50.0),
            phase=params["phase"],
        )
        np.savez(
            OUT / f"{case['name']}.npz",
            plus=np.asarray(h["plus"], dtype=complex),
            cross=np.asarray(h["cross"], dtype=complex),
            spins=np.asarray(spins, dtype=float),  # iota, s1x, s1y, s1z, s2x, s2y, s2z
        )
        print(f"wrote {case['name']}: max|h+| = {np.abs(h['plus']).max():.3e}")


if __name__ == "__main__":
    main()
