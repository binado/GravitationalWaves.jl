"""Generate GWSurrogates reference fixtures from gwsurrogate (NRSur7dq4, NRSur7dq4v2).

Writes ``reference/surrogates/<model>_coorbital.npz``: for each case the mass ratio, the
time-dependent coorbital-frame spins on the ``t_coorb`` grid, the node fit values of one
datapiece, and the coorbital modes returned by gwsurrogate. Also writes
``<model>_dynamics.npz``: dynamics trajectories for several initial conditions and reference
times, plus ODE right-hand sides at probe states and reference times from orbital frequencies.
Finally ``<model>_waveforms.npz`` holds full inertial-frame modes (dimensionless) for several
reference frequencies, output grids and frames, a spin-weighted harmonic table and an end-to-end
strain in SI units.

Run it with ``uv run --project python/surrogates python python/generate_surrogates.py``. With
``--download`` the model files are fetched from Zenodo (CC-BY-4.0, doi:10.5281/zenodo.22257361) and
verified against their SHA-256 first.
"""

import argparse
import hashlib
import urllib.request
from pathlib import Path
from typing import Any

import gwsurrogate.surrogate as gws_top
import numpy as np
from gwsurrogate.new.precessing_surrogate import _utils
from gwtools.harmonics import sYlm
from numpy.typing import NDArray

ROOT = Path(__file__).resolve().parent.parent

# Model name -> datapiece whose node values are dumped to localize fit-evaluation errors.
MODELS: dict[str, str] = {
    "NRSur7dq4": "2_2_Re+",
    "NRSur7dq4v2": "2_2_Re+_sd_0",
}

ZENODO_RECORD = "https://zenodo.org/records/22257361/files"

# Must match MODEL_FILES in packages/GWSurrogates/src/data.jl.
SHA256: dict[str, str] = {
    "NRSur7dq4": "3fc2f58f6ed7a4fff58856b6c01a7dc9155954d2b9fa5d58126a7807826d45d5",
    "NRSur7dq4v2": "62f3afb2baeb837ea19692ff6ca8737808e336cfbd9881c474147ae133efc9de",
}

CASES: list[dict[str, Any]] = [
    {"q": 1.0, "a": (0.0, 0.0), "rate": (0.0, 0.0)},
    {"q": 2.3, "a": (0.5, 0.3), "rate": (9.0, 14.0)},
    {"q": 3.7, "a": (0.78, 0.6), "rate": (25.0, 40.0)},
    {"q": 1.4, "a": (0.2, 0.7), "rate": (3.0, 70.0)},
]


def spin_history(
    t: NDArray[np.float64],
    amplitude: float,
    direction: float,
    rate: float,
    phase: tuple[float, float, float],
) -> NDArray[np.float64]:
    """Build a smooth, norm-bounded spin vector history on the coorbital time grid.

    Parameters
    ----------
    t
        Coorbital time grid, shape ``(N,)``.
    amplitude
        Dimensionless spin magnitude scale.
    direction
        Sign of the offset that tilts the spin towards or away from the orbital axis.
    rate
        Precession-like angular rate over the whole time span.
    phase
        Phase offsets of the three components.

    Returns
    -------
    numpy.ndarray
        Spin vectors with shape ``(N, 3)``.
    """
    s = (t - t[0]) / (t[-1] - t[0])
    chi = np.stack(
        [
            np.cos(rate * s + phase[0]),
            np.sin(rate * s + phase[1]) * (1.0 + 0.3 * s),
            np.cos(0.5 * rate * s + phase[2]) + 0.4 * direction,
        ],
        axis=1,
    )
    chi /= np.linalg.norm(chi, axis=1, keepdims=True)
    return amplitude * (1.0 + 0.1 * np.sin(3 * s))[:, None] * chi


DYNAMICS_CASES: list[dict[str, Any]] = [
    {"q": 1.0, "chiA": (0.0, 0.0, 0.0), "chiB": (0.0, 0.0, 0.0), "t_ref": None},
    {"q": 2.0, "chiA": (0.3, 0.2, 0.4), "chiB": (-0.1, 0.5, -0.3), "t_ref": None},
    {"q": 2.0, "chiA": (0.3, 0.2, 0.4), "chiB": (-0.1, 0.5, -0.3), "t_ref": -100.0},
    {
        "q": 3.2,
        "chiA": (0.5, -0.4, 0.3),
        "chiB": (0.2, 0.3, -0.6),
        "orbphase": 0.7,
        "quat": (0.9, 0.1, -0.3, 0.2),
        "t_ref": 0.0,
    },
    {"q": 1.5, "chiA": (-0.2, 0.1, 0.7), "chiB": (0.4, -0.3, 0.1), "t_ref": -4000.0},
    {"q": 1.5, "chiA": (-0.2, 0.1, 0.7), "chiB": (0.4, -0.3, 0.1), "t_ref": -4260.0},
    {"q": 1.5, "chiA": (-0.2, 0.1, 0.7), "chiB": (0.4, -0.3, 0.1), "t_ref": -4222.0},
    {"q": 1.5, "chiA": (-0.2, 0.1, 0.7), "chiB": (0.4, -0.3, 0.1), "t_ref": -4170.0},
    {"q": 3.9, "chiA": (0.1, 0.6, -0.5), "chiB": (0.0, 0.0, 0.7), "t_ref": 50.0},
]
OMEGA_REFS: list[float] = [0.02, 0.05, 0.1]


def generate_dynamics(name: str, data_dir: Path, out_dir: Path) -> None:
    """Write the dynamics reference fixture of one model.

    Parameters
    ----------
    name
        Model name, a key of ``MODELS``.
    data_dir
        Directory containing ``<name>.h5``.
    out_dir
        Output directory for the ``.npz`` file.
    """
    sur = getattr(gws_top, name)(str(data_dir / f"{name}.h5"))
    dyn = sur._sur_dimless.dynamics_sur
    out: dict[str, Any] = {
        "t_ds": np.asarray(dyn.t),
        "n_cases": np.array(len(DYNAMICS_CASES)),
    }
    for i, case in enumerate(DYNAMICS_CASES):
        quat, orbphase, chi_a, chi_b, _ = dyn(
            case["q"],
            np.array(case["chiA"]),
            np.array(case["chiB"]),
            init_quat=None if "quat" not in case else np.array(case["quat"]),
            init_orbphase=case.get("orbphase", 0.0),
            t_ref=case["t_ref"],
        )
        out[f"case{i}_q"] = np.array(case["q"])
        out[f"case{i}_chiA0"] = np.array(case["chiA"])
        out[f"case{i}_chiB0"] = np.array(case["chiB"])
        out[f"case{i}_orbphase0"] = np.array(case.get("orbphase", 0.0))
        out[f"case{i}_quat0"] = np.array(case.get("quat", (1.0, 0.0, 0.0, 0.0)))
        out[f"case{i}_t_ref"] = np.array(
            np.nan if case["t_ref"] is None else case["t_ref"]
        )
        out[f"case{i}_quat"] = quat
        out[f"case{i}_orbphase"] = orbphase
        out[f"case{i}_chiA"] = chi_a
        out[f"case{i}_chiB"] = chi_b

    # ODE right-hand sides at a probe state, at nodes and at off-node times.
    q = 2.4
    y = np.array([0.8, 0.2, -0.4, 0.4, 1.3, 0.3, -0.2, 0.5, -0.4, 0.1, 0.3])
    y[:4] /= np.linalg.norm(y[:4])
    nodes = [0, 1, 2, 7, 100, 232]
    times = [-4290.0, -4250.0, -3000.5, -100.25, 0.0, 90.0]
    out["probe_q"] = np.array(q)
    out["probe_y"] = y
    out["probe_nodes"] = np.array(nodes)
    out["probe_dydt_nodes"] = np.array(
        [dyn.get_time_deriv_from_index(i, q, y) for i in nodes]
    )
    out["probe_times"] = np.array(times)
    out["probe_dydt_times"] = np.array([dyn.get_time_deriv(t, q, y) for t in times])

    # Reference times from orbital frequencies.
    chi_a, chi_b = np.array([0.3, 0.2, 0.4]), np.array([-0.1, 0.5, -0.3])
    out["omega_refs"] = np.array(OMEGA_REFS)
    out["omega_q"] = np.array(2.0)
    out["omega_chiA"] = chi_a
    out["omega_chiB"] = chi_b
    out["omega_t_refs"] = np.array(
        [dyn._get_t_from_omega(w, 2.0, chi_a, chi_b, 0.0, None) for w in OMEGA_REFS]
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}_dynamics.npz"
    np.savez(path, **out)
    print(f"wrote {path}")


WAVEFORM_CASES: list[dict[str, Any]] = [
    {"q": 2.0, "chiA": (0.3, 0.2, 0.4), "chiB": (-0.1, 0.5, -0.3)},
    {
        "q": 2.0,
        "chiA": (0.3, 0.2, 0.4),
        "chiB": (-0.1, 0.5, -0.3),
        "f_ref": 0.05 / np.pi,
        "dt": 2.0,
        "ell_max": 3,
    },
    {
        "q": 3.1,
        "chiA": (0.5, -0.4, 0.3),
        "chiB": (0.2, 0.3, -0.6),
        "f_ref": 0.05 / np.pi,
        "f_low": 0.04 / np.pi,
        "dt": 1.5,
    },
    {
        "q": 1.7,
        "chiA": (-0.2, 0.1, 0.7),
        "chiB": (0.4, -0.3, 0.1),
        "f_ref": 0.04 / np.pi,
        "times": np.concatenate(
            [np.linspace(-3000.0, -500.0, 200), np.linspace(-499.0, 90.0, 300)]
        ),
    },
    {
        "q": 2.6,
        "chiA": (0.1, 0.6, -0.5),
        "chiB": (0.0, 0.0, 0.7),
        "f_ref": 0.05 / np.pi,
        "dt": 3.0,
        "orbphase": 0.6,
        "quat": (0.9, 0.1, -0.3, 0.2),
    },
    {"q": 1.0, "chiA": (0.0, 0.0, 0.0), "chiB": (0.0, 0.0, 0.0)},
]


def generate_waveforms(name: str, data_dir: Path, out_dir: Path) -> None:
    """Write the full-waveform reference fixture of one model.

    Parameters
    ----------
    name
        Model name, a key of ``MODELS``.
    data_dir
        Directory containing ``<name>.h5``.
    out_dir
        Output directory for the ``.npz`` file.
    """
    sur = getattr(gws_top, name)(str(data_dir / f"{name}.h5"))
    full = sur._sur_dimless
    ell_model = full.coorb_sur.ellMax
    out: dict[str, Any] = {"n_cases": np.array(len(WAVEFORM_CASES))}
    for i, case in enumerate(WAVEFORM_CASES):
        quat = case.get("quat")
        opts = {"init_orbphase": case.get("orbphase", 0.0)}
        if quat is not None:
            opts["init_quat"] = np.array(quat)
        ell_max = case.get("ell_max", ell_model)
        times, h, _ = full(
            [case["q"], np.array(case["chiA"]), np.array(case["chiB"])],
            fM_low=case.get("f_low", 0),
            fM_ref=case.get("f_ref", 0),
            dtM=case.get("dt"),
            timesM=case.get("times"),
            ellMax=ell_max,
            precessing_opts=opts,
        )
        modes = np.array(
            [h[(ell, m)] for ell in range(2, ell_max + 1) for m in range(-ell, ell + 1)]
        )
        out[f"case{i}_q"] = np.array(case["q"])
        out[f"case{i}_chiA0"] = np.array(case["chiA"])
        out[f"case{i}_chiB0"] = np.array(case["chiB"])
        out[f"case{i}_f_ref"] = np.array(case.get("f_ref", 0.0))
        out[f"case{i}_f_low"] = np.array(case.get("f_low", 0.0))
        out[f"case{i}_dt"] = np.array(case.get("dt", np.nan))
        if "times" in case:
            out[f"case{i}_times"] = np.asarray(case["times"])
        out[f"case{i}_ell_max"] = np.array(ell_max)
        out[f"case{i}_orbphase0"] = np.array(case.get("orbphase", 0.0))
        out[f"case{i}_quat0"] = np.array(case.get("quat", (1.0, 0.0, 0.0, 0.0)))
        out[f"case{i}_t"] = np.asarray(times)
        out[f"case{i}_modes"] = modes

    # Spin-weighted spherical harmonics, to pin down the convention of the mode sum.
    thetas = np.array([0.3, 1.1, 2.5])
    phis = np.array([0.0, 0.7, -2.0])
    out["sylm_thetas"] = thetas
    out["sylm_phis"] = phis
    out["sylm"] = np.array(
        [
            [
                sYlm(-2, ell, m, th, ph)
                for ell in range(2, 6)
                for m in range(-ell, ell + 1)
            ]
            for th, ph in zip(thetas, phis)
        ]
    )

    # End to end in SI units: strain at an inclination and reference phase.
    mks = {
        "q": 2.2,
        "chiA": (0.3, 0.2, 0.4),
        "chiB": (-0.1, 0.5, -0.3),
        "M": 60.0,
        "dist_mpc": 400.0,
        "f_ref": 25.0,
        "dt": 1.0 / 4096,
        "inclination": 0.7,
        "phi_ref": 0.3,
    }
    domain, strain, _ = sur(
        mks["q"],
        np.array(mks["chiA"]),
        np.array(mks["chiB"]),
        M=mks["M"],
        dist_mpc=mks["dist_mpc"],
        f_low=0,
        f_ref=mks["f_ref"],
        dt=mks["dt"],
        inclination=mks["inclination"],
        phi_ref=mks["phi_ref"],
        units="mks",
    )
    for key, value in mks.items():
        out[f"mks_{key}"] = np.array(value)
    out["mks_t"] = np.asarray(domain)
    out["mks_h"] = np.asarray(strain)

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}_waveforms.npz"
    np.savez(path, **out)
    print(f"wrote {path}")


def download(name: str, data_dir: Path) -> None:
    """Fetch the HDF5 file of one model from Zenodo and verify its checksum.

    Parameters
    ----------
    name
        Model name, a key of ``SHA256``.
    data_dir
        Directory the file is written to.

    Raises
    ------
    ValueError
        If the SHA-256 of the downloaded file does not match.
    """
    path = data_dir / f"{name}.h5"
    if path.exists():
        return
    data_dir.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(".part")
    urllib.request.urlretrieve(f"{ZENODO_RECORD}/{name}.h5", part)
    checksum = hashlib.sha256(part.read_bytes()).hexdigest()
    if checksum != SHA256[name]:
        part.unlink()
        raise ValueError(
            f"checksum mismatch for {name}.h5: expected {SHA256[name]}, got {checksum}"
        )
    part.rename(path)


def generate(name: str, data_dir: Path, out_dir: Path) -> None:
    """Write the coorbital reference fixture of one model.

    Parameters
    ----------
    name
        Model name, a key of ``MODELS``.
    data_dir
        Directory containing ``<name>.h5``.
    out_dir
        Output directory for the ``.npz`` file.
    """
    probe = MODELS[name]
    sur = getattr(gws_top, name)(str(data_dir / f"{name}.h5"))
    coorb = sur._sur_dimless.coorb_sur
    ell_max = coorb.ellMax
    t = np.asarray(coorb.t)
    out: dict[str, Any] = {
        "t_coorb": t,
        "ell_max": np.array(ell_max),
        "n_cases": np.array(len(CASES)),
    }
    for i, case in enumerate(CASES):
        q = case["q"]
        chi_a = spin_history(t, case["a"][0], 1.0, case["rate"][0], (0.1, 0.7, 1.3))
        chi_b = spin_history(t, case["a"][1], -1.0, case["rate"][1], (0.9, 0.2, 2.1))
        mid = len(t) // 2
        fit_params = coorb._get_fit_params(np.array([q, *chi_a[mid], *chi_b[mid]]))
        out[f"case{i}_q"] = np.array(q)
        out[f"case{i}_chiA"] = chi_a
        out[f"case{i}_chiB"] = chi_b
        out[f"case{i}_modes"] = coorb(q, chi_a, chi_b, ellMax=ell_max)
        out[f"case{i}_fit_params_mid"] = fit_params

    # Node values of one datapiece at case 1, to localize errors in the fit evaluation.
    data = coorb.data[probe]
    q = CASES[1]["q"]
    chi_a, chi_b = out["case1_chiA"], out["case1_chiB"]
    nodes = []
    for orders, coefs, ni in zip(data["orders"], data["coefs"], data["nodeIndices"]):
        x = np.empty(7)
        x[0] = q
        x[1:4] = chi_a[ni]
        x[4:7] = chi_b[ni]
        nodes.append(
            _utils.eval_fit(
                orders, coefs, coorb._get_fit_params(x), *coorb._fit_settings
            )
        )
    out["probe_nodes"] = np.array(nodes)

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}_coorbital.npz"
    np.savez(path, **out)
    print(f"wrote {path}")


def main() -> None:
    """Parse arguments and generate the fixtures of every available model."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path, required=True, help="directory with <model>.h5"
    )
    parser.add_argument("--out", type=Path, default=ROOT / "reference" / "surrogates")
    parser.add_argument(
        "--download", action="store_true", help="fetch missing model files from Zenodo"
    )
    args = parser.parse_args()
    for name in MODELS:
        if args.download:
            download(name, args.data_dir)
        if (args.data_dir / f"{name}.h5").exists():
            generate(name, args.data_dir, args.out)
            generate_dynamics(name, args.data_dir, args.out)
            generate_waveforms(name, args.data_dir, args.out)
        else:
            print(f"skipping {name}: {args.data_dir / (name + '.h5')} not found")


if __name__ == "__main__":
    main()
