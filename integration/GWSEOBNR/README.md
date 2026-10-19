# SEOBNR acceptance

This environment is separate from the workspace and default test targets. It
checks all nine specified physical cases on `0:0.0625:2048` Hz, raw and optimized
mismatch, component matching, directional ForwardDiff derivatives, five-point
finite differences at five perturbations, extrinsic scaling, tighter solver
convergence, and an offset H1/L1 fixed-PSD likelihood gradient. Intrinsic finite
differences use the internal signed-spin continuation just outside an envelope
boundary; the public API continues to reject out-of-envelope draws.

```
just fixtures-seobnr
just test-integration-seobnr
```

The Python reference environment requires GSL, SWIG, and a C compiler. Consult
`python/seobnr/README.md` for the pinned pygsl-lite build compatibility flag.
Normal package tests use compact committed snapshots and need no Python or LAL.

The acceptance runner intentionally returns a failing status for unmet gates.
It saves partial results after each case and final results to
`reference/seobnr/diagnostics.toml`, including solver steps, branch choices,
reference and attachment times, NQC coefficients, spectral errors, derivative
errors, convergence, warmed timings, and allocations. The analytic positive PSD
is proportional to `1 + (40/f)^4 + (f/200)^2` above 20 Hz, with finite positive
excluded bins. Diagnostic time/phase optimization never changes primary gates.

**Acceptance is pending.** See `packages/GWSEOBNR/ACCEPTANCE.md` and its committed
measurement record. NQC matching and several spectral L² comparisons remain
outside the original tolerances. These failures are not skipped or relaxed.

Exploratory native Enzyme reverse-mode results and reproduction instructions are
in [ENZYME.md](ENZYME.md). These probes do not establish reverse-mode support.
