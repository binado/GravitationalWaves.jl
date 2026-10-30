# GWSeries

MIT-licensed `TimeSeries` matrices for uniform time segments and `FrequencySeries`
matrices for complete, zero-origin, one-sided spectra, built on
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl). Rows are samples
or frequency bins (`GPSTime` or `Freq` dimensions); columns are ordered channels (`Chan`).
Ordinary matrices, vectors (reshaped to one column), and views can be wrapped without
copying samples. Metadata precision is independent of sample precision.

```julia
using GWSeries
x = TimeSeries(randn(8192, 2); epoch=1126259460.125,
               sampling_frequency=4096, channels=(:H1, :L1))
array, meta = flatten(x)
y = TimeSeries(2 .* array; meta...)
parent(y)                       # ordinary matrix
grid_metadata(y)                # canonical constructor keywords
times(y)                        # Base range of absolute GPS sample times
```

`grid_metadata(TimeSeries)` contains `epoch`, `sampling_interval`, and `channels`.
`grid_metadata(FrequencySeries)` contains `epoch`, `frequency_spacing`, `channels`, and
`time_samples`. (DimensionalData's own `metadata` is unused and returns `NoMetadata()`.)
The grid lives in the dimension lookups: `lookup(x, GPSTime)` wraps a `TimeGrid` of
absolute sample times and `lookup(h, Freq)` a `SpectralGrid` from DC, which also records
the epoch and original time length. `times(x)` and `frequencies(h)` return their raw Base
ranges (`range(epoch; step, length)` and `range(0; step = df, length)`), so code that
computes with coordinates gets ordinary Base range arithmetic. Sample count is derived from storage,
so reconstruction with different row counts changes duration/extent. Constructors
validate intrinsic metadata, one-based storage axes, and unique channel names.

## Array operations and selection

Series are `AbstractDimArray`s, so DimensionalData selectors, broadcasting, reductions,
plotting recipes, and Tables interop apply. A result is still a series only when it
keeps a valid grid; otherwise it is a plain `DimArray` that keeps the labels but not
the series guarantees.

| Operation | Result |
|---|---|
| same-grid `+`, `-`, scalar `*` and `/`, broadcasting, `map`, `copy`, `similar`, `cumsum(; dims=1)`, contiguous or strided time slices, channel selection | series |
| `reverse`, `circshift`, `diff`, partial frequency bands, reductions (`sum(; dims=1)`), `x * A`, `x * x'`, `permutedims`, repeated channels, `cumsum(; dims=2)` | `DimArray` |
| grid, epoch, or channel-order mismatch between operands | throws |

```julia
2 .* x                          # TimeSeries, unchanged grid
x[Chan = At(:H1)]               # DimArray vector with GPS times
x[Chan = At([:L1, :H1])]        # TimeSeries, channels in the requested order
x[101:200, :]                   # copy, advanced epoch
x[GPSTime = t0 .. t0 + 1]       # closed interval of absolute GPS times
x[1:2:end, :]                   # TimeSeries with doubled sampling interval
view(x, 101:200, :)             # TimeSeries wrapping shared storage
h[Freq = 20 .. 1024]            # DimArray: a band is not a complete spectrum
```

`At`, `Near`, and `..` are re-exported; load DimensionalData for its other selectors.
Positional and keyword indices cannot be mixed in one call (`x[GPSTime = 1:4, Chan = At(:H1)]`
works). Reversed or irregular time selections return `DimArray`s. Frequency selections
preserve `FrequencySeries` only for the full axis; partial bands drop the epoch and time
length. In-place reordering (`reverse!`, `sort!`) cannot change the type, so avoid it on
series. `crop(x; start, stop)` and `pad(x; before, after)` apply only to TimeSeries. Crop
uses grid-aligned half-open boundaries; padding before shifts the epoch.

Operands are checked before executing on parent storage: broadcasting (including into a
destination with `.=`), `+`, `-`, `map`, `copyto!`, and assignment of a series
(`out[:, :] = y`) compare grids. Use `flatten` around repeated numerical work to avoid
repeated checks. Static matrices combine with `+`/`-` into unlabelled static arrays;
broadcasting them with any DimensionalData array is unsupported upstream.

`frequency_bins(h; fmin, fmax)` returns indices of existing bins inside inclusive bounds,
clipped to the available grid, without rounding bounds onto a new grid. Disjoint bounds
produce an empty range; nonfinite or reversed bounds throw. Keep these indices separately
when selecting a band.

Coordinate comparison retains exact fast paths and a few-ulp tolerance: four ulps
of absolute coordinates plus eight ulps of displacement and spacing. Ambiguous
comparisons whose tolerance reaches half a sample/bin are rejected. BigFloat
metadata is available for GPS epochs requiring sub-microsecond precision.
DimensionalData itself compares numeric lookups with a relative `≈`, which would accept
GPS epochs seconds apart; GWSeries replaces that comparison for `GPSTime` and `Freq`
grids through a DimensionalData-internal hook (hence the narrow `DimensionalData` compat
bound), so broadcasting with `strict_broadcast` (the default) applies the tolerance above.

## Physical AbstractFFTs transforms

GWSeries uses the AbstractFFTs API and does not load an FFT backend. Load one
supporting your parent storage (for example FFTW for CPU arrays):

```julia
using AbstractFFTs, FFTW
h = rfft(x)                    # FrequencySeries: dt*rfft(parent(x), 1)
z = irfft(h)                   # TimeSeries, original epoch, channels and length
z = irfft(h, nsamples(x))      # explicit length, checked against the recorded one
h = to_frequencydomain(x)      # real samples only
z = to_timedomain(h)
```

GWSeries depends only on the AbstractFFTs interface. If no backend is loaded for the
parent storage, series transforms and plans throw an `ArgumentError` naming the fix
(`using FFTW` for CPU arrays).

Series transforms act only on dimension one, preserve epoch and channels, and
use physical normalization. Transforms of `parent(x)` retain standard backend
conventions. Real transforms support odd and even lengths; the original length is always required
for an inverse, including inverse plans. A spectrum must have `N÷2+1` bins.
Complex/two-sided transforms operate on `parent(series)` without series metadata.
No absolute-GPS phase factor is applied. Imaginary DC/even-Nyquist components
follow the backend's real-inverse convention.

```julia
forward = plan_to_frequencydomain(x)
h = transform(forward, x)
inverse = plan_to_timedomain(h)
grid_metadata(forward)          # output constructor keywords
transform!(h, forward, x)
raw_h = transform(forward, parent(x))
```

Prepared transforms use AbstractFFTs plans, physical scaling, and reusable private
scratch storage. Inputs are preserved. Plans are not reentrant; use separate
plans for concurrent tasks. Raw execution accepts prevalidated arrays; wrapped
execution also checks metadata. DimensionalData's own AbstractFFTs methods (`plan_rfft`,
`fft`, ...) are rejected for series because they apply an absolute-epoch phase and drop the
grid; use the physical methods above. Backend-specific preparation keywords pass to
AbstractFFTs. Unsupported storage never silently moves to the CPU. Type-changing
in-place series FFTs are rejected; use a plan with reusable output storage.

## Fixed windows and compiled work

`PlanckWindow(x; start_width, end_width)` prepares a fixed taper. Both widths are
required, in seconds. Alternatively prepare with `PlanckWindow(; nsamples,
sampling_interval, start_width, end_width)`. Windows ignore epoch/channel labels
and retain the existing continuous Planck transition and floating-point precision.

`apply_window(window, x)` preserves series metadata; `apply_window!` supports
series destinations. Raw array methods remain unchecked and support higher-dimensional
arrays. Prepare windows outside differentiation.

AD acceptance covers fixed metadata and numerical samples, including the
`flatten → raw kernel → reconstruction` boundary. FFT differentiation is outside
this contract. ArrayInterface traits come from DimensionalData's extension and forward to
the parent storage.

## Inverting a selected band

A band is a `DimArray`, not an independently invertible spectrum. Restore its bins
explicitly to produce a band-limited signal; zeros do not recover discarded values:

```julia
bins = frequency_bins(h; fmin=20.1, fmax=1024.0)
band = view(h, bins, :)
full = copy(h)
fill!(full, 0)
full[bins, :] = band
t_bandlimited = irfft(full)
```

A spectrum from `rfft` records its original length (`time_samples(h)`), which removes the
even/odd ambiguity of the inverse: lengths `N` and `N + 1` can share a bin count, and
`irfft(h, N + 1)` throws for a spectrum recorded at `N`. A `FrequencySeries` built from
elsewhere, such as a waveform, has `time_samples(h) === nothing` and needs the length
passed explicitly. Changing the number of bins while holding Δf fixed changes the time
sampling implied by a chosen inverse length.

### Differentiating the forward Fourier transform

The allocating `to_frequencydomain(::TimeSeries)` / `rfft(::TimeSeries)` path
supports sample derivatives with ForwardDiff, or Mooncake together with
ChainRulesCore. Load an FFT backend (for example FFTW) separately. Metadata,
window construction, and transform configuration remain fixed. ForwardDiff
transforms the primal and tangent arrays with the same FFT; the Mooncake extension
reuses AbstractFFTs' ChainRules pullback, including one-sided endpoint weights.
Mooncake's bridge currently supports dense Float32/Float64 sample matrices with
matching spacing precision. Inverse transforms and mutable prepared-plan AD are
not covered by this contract. See `just test-integration-series` for validation.
