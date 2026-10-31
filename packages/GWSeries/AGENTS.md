# GWSeries agent guidance

GWSeries provides DimensionalData-backed matrices for uniformly sampled time series and complete one-sided spectra.

- Keep storage in `(sample, channel)` order; `flatten` must return parent storage
  and reconstructible `grid_metadata` keywords without losing metadata precision.
- Grid metadata lives only in the lookups (`TimeGrid`, `SpectralGrid`,
  `ChannelNames`); series structs hold just `data` and `dims`. `DD.rebuild` on a
  series is the single gate deciding whether a result stays a series: it must
  return a series only for canonical grid dims with unique channels, and a plain
  `DimArray` otherwise. Grid `getindex` methods decide which slices keep a grid
  (contiguous or strided time slices do; partial frequency bands do not).
- Grids wrap the Base range they represent and delegate all values to it.
  `times`/`frequencies` must return that raw range, never the grid type: Base
  range arithmetic on custom `AbstractRange`s rebuilds via `range` and can shift
  coordinates by an ulp, changing downstream kernel results.
- Grid coordinates are any `Real` (epochs may be dual numbers under ForwardDiff).
  Keep grid fields concrete: a `Nothing`-typed field nested in the dims crashed
  Mooncake's rule compilation, hence the zero sentinel for unrecorded lengths.
- Add explicit methods for DD operations that keep dims while reordering or
  mixing samples (`circshift`, `diff`, `cumsum` across channels), and for paths
  that skip lookup-value checks (`map`, `copyto!`, series `setindex!`).
- Keep private DD 0.30 dependencies isolated in `src/dd_compat.jl`:
  `Dimensions._comparedims2` is required because tuple comparison bypasses public
  pairwise dispatch; epoch-shift, predicate, and singleton-broadcast regressions
  guard its tolerance and behavior. `DD.StandardIndices` routes terminal assignment
  after selector conversion; `DD._DimIndicesAmb` resolves overlapping selector and
  empty-index-array methods. Mixed labelled assignment and empty-index regressions,
  plus Aqua, cover these aliases. Keep specializations on GWSeries-owned types,
  retain the narrow compatibility bound, and recheck all adapters on DD upgrades.
  Use public `DD.Dimensions.comparedims` for ordinary operand validation.
- Methods overriding DD's AbstractFFTs extension must mirror its fixed arities;
  `Vararg` rejection methods are ambiguous with it. Run Aqua after FFT changes.
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
