# SEOBNR oracle (separate from the bilby environment)

Pinned pySEOBNR 0.3.7 and qnm 0.4.4, with transitive versions and archive hashes
in `uv.lock`. Requires a C compiler, SWIG, and GSL development libraries.

```
uv sync --project python/seobnr --locked
uv run --project python/seobnr python python/seobnr/generate.py --components
uv run --project python/seobnr python python/seobnr/export_qnm.py
uv run --project python/seobnr python python/seobnr/generate.py --waves
```

The pinned pygsl-lite 0.1.8 source may require
`CFLAGS=-DPyInt_FromLong=PyLong_FromLong` when built with SWIG 4.5.1 on macOS.
This Python 3 C API compatibility flag does not change the physical equations.
The setuptools build constraint matches pySEOBNR's shipped constraint.

Post-adiabatic evolution is disabled and only (2,2) output requested. This does
not reduce the multipoles used by the radiation-reaction flux. Waveforms and
trajectories are written under gitignored `reference/seobnr`. Small component
snapshots and exact QNM interpolation coefficients are committed with the package.

After generating large mode fixtures, apply the same continuous attachment-peak
convention and physical conditioning, and regenerate the compact short case:

```
uv run --project python/seobnr python python/seobnr/condition.py
uv run --project python/seobnr python python/seobnr/short_snapshot.py
```

The raw pySEOBNR output remains available in the original `.npz` files. The
conditioner recovers the continuous (2,2) attachment peak from its zero-slope
matching target, undoing the sampled upstream output-time shift. It does not
optimize time or phase against Julia. QNM sequence hashes are recorded in
`packages/GWSEOBNR/qnm-provenance.toml` alongside the bundled interpolants.

To reproduce the scalar equation ports, extract the archive identified in
`packages/GWSEOBNR/provenance.toml`, then run from the repository root:

```
python3 python/seobnr/port_equations.py /path/to/pyseobnr-0.3.7
runic --inplace packages/GWSEOBNR/src
```

The generator verifies committed source-file checksums before writing outputs.
