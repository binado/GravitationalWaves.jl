# GWParameterConversion agent guidance

GWParameterConversion provides dependency-free, AD-compatible mass, spin, and frame conversions following bilby.

- Preserve bilby names, units, and floating-point operation order, including
  `mass_ratio = mass_2 / mass_1 ≤ 1` and `Float64(π)` when mirroring NumPy comparisons.
- Keep numerical arguments generic for real-valued physics, preserving dual numbers
  and BigFloat precision. Avoid introducing runtime dependencies for conversions.
- Keep target keys, consistency flags, and tolerances fixed during differentiation.
  `convert_parameters` retains supplied keys; disabling consistency checks requires
  consistent redundant inputs. Preserve documented physical-domain preconditions.
- Preserve finite equal-mass derivatives in coordinate transforms and their log
  Jacobians. Test boundary values separately from derivatives at nonsmooth spin
  switches, zero norms, undefined azimuths, or singular inverse mass coordinates.
- Generate bilby fixtures from the repository root with `uv sync --project python`
  and `uv run --project python python python/generate.py`. Missing fixtures warn
  locally and fail under `CI=true`; the package suite includes ForwardDiff checks.
- Keep optional Jacobian/inverse integrations independent; add no Jacobian method
  for `convert_parameters`, which retains inputs and adds derived fields.

See the [README](README.md) for domain, consistency, and coordinate-change contracts,
and the [repository README](../../README.md) for reference-fixture setup.
