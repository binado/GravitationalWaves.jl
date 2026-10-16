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
    (OUT / "conversion.json").write_text(json.dumps({"masses": cases}, indent=1))
    print(f"wrote conversion.json: {len(cases)} mass cases")


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
