# GWNoise agent guidance

GWNoise reads one-sided PSD/ASD tables and samples finite stationary Gaussian noise segments.

- Keep numerical APIs on plain arrays without a GWSeries dependency: data have
  shape `(nf, nc)`, PSDs are `(nf, nc)` real or `(nc, nc, nf)` Hermitian.
- Use full-length PSD arrays with zeros marking inactive bins, never sliced bands.
  Physical Fourier coefficients satisfy `E[n n†] = S / (2Δf)`.
- Preserve `band_psd` endpoint behavior: zero DC and the last bin, assuming an
  even-length transform. An odd-length transform also loses its last real bin.
- Return `W = L⁻¹`, where `L L† = S`, from `whitening_factor`, with zero inactive
  bins. Invalid active blocks (non-finite, non-Hermitian, or non-positive-definite)
  must throw `ArgumentError` during preparation.
- Keep table units explicit (Hz, strain²/Hz or ASD strain/√Hz); square ASD before
  interpolation and reject out-of-coverage evaluation unless zero fill is requested.
- Draw random numbers for every bin so streams stay aligned across bands; never
  reseed caller RNGs. PSD estimation and continuous streaming remain out of scope.

See the [README](README.md) for table formats, interpolation, and noise normalization.
