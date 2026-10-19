"""Capture Ripple's XAS phase systems using NumPy (no JAX dependency).

Usage: python3 generate_fixtures.py /path/to/pinned/ripple/waveforms
The directory must contain the two files from RIPPLE_REVISION below.
Only definitions/constants are loaded; this does not run Ripple's waveform driver.
"""

import ast
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

RIPPLE_REVISION = "58c53024718697e0430b666ef5dd40ff3c2e029d"
SOURCE_HASHES = {
    "IMRPhenomXAS.py": "303bf59c3251aa108e57ed9b7a4f829fa8ff7298b4ce18adbf4b9a2f7d20c0b3",
    "IMRPhenomX_utils.py": "3eed503072d3a51eef7fce5e357d999637b99de99c9fce397288c1c6d51d2f7d",
}


class CaptureNumpy:
    def __init__(self):
        self.system = None
        self.linalg = SimpleNamespace(solve=self.solve)

    def __getattr__(self, name):
        return getattr(np, name)

    def solve(self, a, b):
        self.system = (a.copy(), b.copy())
        return np.linalg.solve(a, b)


def definitions(path, namespace):
    source = path.read_text()
    digest = hashlib.sha256(source.encode()).hexdigest()
    if digest != SOURCE_HASHES[path.name]:
        raise ValueError(f"{path.name} does not match pinned revision {RIPPLE_REVISION}")
    tree = ast.parse(source)
    tree.body = [
        node for node in tree.body
        if isinstance(node, (ast.Assign, ast.FunctionDef))
    ]
    exec(compile(tree, str(path), "exec"), namespace)
    return digest


def main():
    source = Path(sys.argv[1])
    capture = CaptureNumpy()
    # gt cancels from the dimensionless frequencies and mass fractions here.
    # Set it to one to keep these fixtures independent of unit constants.
    common = dict(jnp=capture, gt=1.0, PI=np.pi, Array=np.ndarray, Tuple=tuple)
    utils = common.copy()
    utility_hash = definitions(source / "IMRPhenomX_utils.py", utils)
    model = dict(common, EulerGamma=np.euler_gamma,
                 IMRPhenomX_utils=SimpleNamespace(**utils))
    model_hash = definitions(source / "IMRPhenomXAS.py", model)
    coeffs = utils["PhenomX_phase_coeff_table"]
    lines = [
        f'ripple_revision = "{RIPPLE_REVISION}"',
        f'model_sha256 = "{model_hash}"',
        f'utilities_sha256 = "{utility_hash}"',
        'oracle = "Pinned Ripple formulas evaluated with NumPy float64"',
    ]
    draws = [
        ("unequal", [36.0, 29.0, 0.3, -0.2]),
        ("equal", [30.0, 30.0, 0.0, 0.0]),
        ("equal_unequal_spin", [30.0, 30.0, 0.6, -0.3]),
        ("aligned", [50.0, 10.0, 0.9, 0.7]),
        ("antialigned", [50.0, 10.0, -0.8, -0.6]),
        ("high_ratio", [100.0, 1.0, 0.2, -0.1]),
    ]
    for name, theta in draws:
        theta = np.array(theta)
        rd, damp, meco, isco = utils["get_cutoff_fMs"](*theta)
        f2base = 0.6 * (0.5 * rd + isco)
        f1 = meco - (f2base - meco) * 0.03
        model["get_inspiral_phase"](np.float64(0.02), theta, coeffs)
        inspiral = capture.system
        # Complex step is used only to supply the intermediate endpoint value;
        # the Julia prototype treats every RHS entry as an independent input.
        din = model["get_inspiral_phase"](f1 + 1j * 1e-20, theta, coeffs).imag / 1e-20
        _, (cl, rd0) = model["get_mergerringdown_raw_phase"](
            np.float64(0.08), theta, coeffs
        )
        ringdown = capture.system
        model["get_intermediate_raw_phase"](
            np.float64(0.04), theta, coeffs, din, rd0, cl
        )
        intermediate = capture.system
        for region, (a, b) in zip(
            ("inspiral", "intermediate", "ringdown"),
            (inspiral, intermediate, ringdown),
        ):
            assert np.isfinite(a).all() and np.isfinite(b).all()
            lines.extend([
                "", "[[case]]", f'name = "{name}/{region}"',
                f'region = "{region}"',
                "theta = " + json.dumps(theta.tolist()),
                "frequencies = " + json.dumps([meco, isco, rd, damp]),
                "rhs = " + json.dumps(b.tolist()),
                "matrix = " + json.dumps(a.tolist()),
                "solution = " + json.dumps(np.linalg.solve(a, b).tolist()),
            ])
    target = Path(__file__).with_name("fixtures.toml")
    target.write_text("\n".join(lines) + "\n")
    print(f"Wrote {len(draws) * 3} systems to {target}")


if __name__ == "__main__":
    main()
