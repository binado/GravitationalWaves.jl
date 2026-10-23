# GWDetectors agent guidance

GWDetectors provides differentiable interferometer geometry, LAL-convention timing and antenna factors, and FD projection.

- Preserve pinned LAL Cartesian constants, WGS-84 geometry, and antenna conventions.
  Do not renormalize Cartesian arms or force Float64; use StaticArrays and preserve
  promoted numeric types. Arm lengths are metadata for the long-wavelength response.
- Keep delay signs (detector minus geocenter is `-dot(vertex, direction)/c`) and
  LAL J2000 GMST conventions. Preserve whole/fractional GPS separation and the
  pinned leap-table behavior documented in the README.
- Projection consumes plus/cross series and returns detector series on uniform
  integer-indexed frequency bins. Preserve LAL/bilby phase operation order; subtract
  destination GPS epoch before adding the source offset to retain precision.
- Keep metadata and frequency/epoch configuration fixed for AD; raw projection
  accepts prevalidated arrays. Geometry validation may be bypassed with `check=false`
  for prevalidated differentiated inputs; retain size/type checks.
- Generate fixtures from the repository root with `uv sync --project python` and
  `uv run --project python python python/generate_detectors.py`. Missing fixtures
  warn locally and fail under `CI=true`.
- Run `just test-integration-detectors` for ForwardDiff/Mooncake acceptance; test
  GPS derivatives away from rounding boundaries and leap transitions.

See the [README](README.md) for projection precision, supported timing domains, and geometry conventions.
