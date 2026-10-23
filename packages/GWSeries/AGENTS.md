# GWSeries agent guidance

GWSeries provides matrix wrappers for uniformly sampled time series and complete one-sided spectra.

- Keep storage in `(sample, channel)` order; `flatten` must return parent storage
  and reconstructible keyword metadata without losing metadata precision.
- Preserve metadata for same-shape `+`/`-`, scalar scaling, and broadcasting.
  Matrix products, `similar`, and sample-reordering operations return plain arrays.
  Contiguous time slices update epoch and channel order; partial frequency slices
  return plain arrays, since they are not complete spectra.
- FrequencySeries starts at DC and records the original `time_samples` after
  `rfft`; otherwise that value is `nothing`. Inverse transforms use the recorded
  length or require it explicitly, preserving odd/even distinctions.
- Use AbstractFFTs with a separately loaded backend. Physical transforms act on
  dimension one as `dt * rfft`, with no absolute-epoch phase; preserve epoch and
  channels. Never silently transfer unsupported storage to the CPU.
- Keep metadata, windows, and transform configuration fixed during AD. Preserve
  the raw-array/reconstruction boundary; prepared plans are not reentrant.
  Forward sample AD is supported; inverse and mutable prepared-plan AD are outside
  the acceptance contract.
- Run `just test-integration-series` for Fourier or AD changes and test downstream
  waveform, interface, detector, and shared-test consumers for wrapper changes.

See the [README](README.md) for slicing, normalization, plan, and AD contracts.
