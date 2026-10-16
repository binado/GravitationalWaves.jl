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


def main() -> None:
    spec = json.loads((HERE / "cases.json").read_text())
    df = spec["frequency"]["delta_f"]
    f_max = spec["frequency"]["f_max"]
    frequency = np.arange(int(round(f_max / df)) + 1) * df

    OUT.mkdir(exist_ok=True)
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
