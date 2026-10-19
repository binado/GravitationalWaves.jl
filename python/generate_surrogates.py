"""Generate GWSurrogates reference fixtures from gwsurrogate (NRSur7dq4, NRSur7dq4v2).

Writes ``reference/surrogates/<model>_coorbital.npz``: for each case the mass ratio, the
time-dependent coorbital-frame spins on the ``t_coorb`` grid, the node fit values of one
datapiece, and the coorbital modes returned by gwsurrogate.

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
        else:
            print(f"skipping {name}: {args.data_dir / (name + '.h5')} not found")


if __name__ == "__main__":
    main()
