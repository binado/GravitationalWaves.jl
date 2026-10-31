# GWMatchedFiltering agent guidance

GWMatchedFiltering computes whitened frequency-domain statistics on plain arrays without runtime dependencies.

- Use `(nf, nc)` data (or vectors for one channel) and `W = L⁻¹`, `L L† = S`:
  `(nf, nc)` real diagonal weights or `(nc, nc, nf)` lower-triangular blocks.
  Zero weights mark inactive bins; conjugate the first argument in
  `⟨a,b⟩ = 4Δf Σₖ (W aₖ)†(W bₖ)`.
- Keep structure dispatch in `whiten` through `filter_structure`; derive every
  statistic from `whitened_contraction` rather than duplicating diagonal/block logic.
- Kernels assume prevalidated arrays with identical shapes for `a` and `b`.
  Matched-filter SNR requires a finite nonzero template: call `check_template`
  explicitly on the host.
- Keep `W` and `df` fixed for AD; preserve the generic array-expression fallback.
  Any future host fast path must require **all** of `a`, `b`, and `W` to be host
  arrays, including when weights are constant and strain arrays are traced.
- Use `just test-integration-inference` for ForwardDiff/Mooncake and composed
  waveform/projection/likelihood acceptance; use `just benchmark-inference` for
  performance changes and storage-layout comparisons.

See the [README](README.md), [testing guide](../../README.md#development-and-tests), and
[benchmark guide](../../benchmarks/GWMatchedFiltering/README.md).
