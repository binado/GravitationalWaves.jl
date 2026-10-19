"""Export exact piecewise-cubic qnm 0.4.4 interpolants used by pySEOBNR."""
import hashlib
from pathlib import Path
import qnm
import scipy.interpolate
from generate import verify_versions, ROOT

verify_versions()
lines = ["# SPDX-License-Identifier: MIT", "# qnm 0.4.4 (Leo C. Stein); see NOTICE.md and qnm-LICENSE.",
         "# Exact piecewise cubic coefficients, exported by python/seobnr/export_qnm.py."]
snapshots = []
for m, name in ((2, "POSITIVE"), (-2, "NEGATIVE")):
    seq = qnm.modes_cache(s=-2, l=2, m=m, n=0)
    for part, interp in (("REAL", seq._interp_o_r), ("IMAG", seq._interp_o_i)):
        poly = scipy.interpolate.PPoly.from_spline(interp._eval_args)
        # Repeated knots have zero width; retain only genuine intervals.
        mask = poly.x[1:] > poly.x[:-1]
        knots = [float(x) for x in poly.x[:-1][mask]] + [float(poly.x[-1])]
        coeffs = poly.c[:, mask]
        lines += [f"const QNM_{name}_{part}_KNOTS = {knots!r}",
                  f"const QNM_{name}_{part}_COEFFS = ["]
        lines += ["    " + " ".join(repr(float(x)) for x in row) + ";" for row in coeffs.T]
        lines += ["]"]
    for a in (0., .1, .4, .7, .9, .98):
        omega = seq(a=a, interp_only=True)[0]
        snapshots.append(dict(m=m, spin=a, real=float(omega.real), imag=float(omega.imag)))
path = ROOT / "packages/GWSEOBNR/src/qnm_data.jl"
path.write_text("\n".join(lines) + "\n")
import tomlkit
(ROOT / "packages/GWSEOBNR/test/qnm_oracle.toml").write_text(tomlkit.dumps(dict(cases=snapshots)))
print(path, hashlib.sha256(path.read_bytes()).hexdigest(), flush=True)
