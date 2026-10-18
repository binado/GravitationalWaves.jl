"""Generate scalar bilby and NumPy correlated-metric references (no LAL calls)."""

import json
from pathlib import Path

import bilby
import numpy as np
from bilby.gw.utils import matched_filter_snr, noise_weighted_inner_product


def main():
    f = np.arange(5) * 0.5
    h = np.array([[0, 0], [1 + 2j, 2 - 1j], [2j, -1], [-1, 1j], [0.5, 2]])
    d = np.array([[0, 0], [2 - 1j, -1 + 3j], [1j, 2], [1 + 2j, -1j], [1, -0.5j]])
    psds = np.array([[0, 0], [2, 3], [3, 4], [4, 5], [5, 6]])
    band = slice(1, 5)
    duration = 1 / (f[1] - f[0])
    overlap = sum(
        noise_weighted_inner_product(h[band, i], d[band, i], psds[band, i], duration)
        for i in range(2)
    )
    norm = sum(
        noise_weighted_inner_product(h[band, i], h[band, i], psds[band, i], duration).real
        for i in range(2)
    )
    single = matched_filter_snr(h[band, 0], d[band, 0], psds[band, 0], duration)
    spectra = np.zeros((5, 2, 2), dtype=complex)
    spectra[:, 0, 0] = psds[:, 0]
    spectra[:, 1, 1] = psds[:, 1]
    spectra[band, 0, 1] = 0.2 + 0.3j
    spectra[band, 1, 0] = 0.2 - 0.3j
    correlated = 4 / duration * sum(
        np.vdot(h[k], np.linalg.solve(spectra[k], d[k])) for k in range(1, 5)
    )
    result = {
        "bilby_version": bilby.__version__,
        "numpy_version": np.__version__,
        "diagonal_overlap": [float(overlap.real), float(overlap.imag)],
        "diagonal_norm": float(norm),
        "single_matched": [float(single.real), float(single.imag)],
        "correlated_overlap": [float(correlated.real), float(correlated.imag)],
    }
    destination = Path(__file__).resolve().parents[1] / "reference" / "inference.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
