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
