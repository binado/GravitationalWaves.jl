# GWSeries

MIT-licensed `TimeSeries` matrices for uniform time segments and `FrequencySeries`
matrices for complete, zero-origin, one-sided spectra. Rows
are samples or frequency bins; columns are ordered channels. Ordinary matrices,
vectors (reshaped to one column), and views can be wrapped without copying samples.
Metadata precision is independent of sample precision.

```julia
using GWSeries
x = TimeSeries(randn(8192, 2); epoch=1126259460.125,
               sampling_frequency=4096, channels=(:H1, :L1))
array, meta = flatten(x)
y = TimeSeries(2 .* array; meta...)
parent(y)                       # ordinary matrix
metadata(y)                     # canonical constructor keywords
```

`metadata(TimeSeries)` contains `epoch`, `sampling_interval`, and `channels`.
`metadata(FrequencySeries)` contains only `epoch`, `frequency_spacing`, and `channels`.
Frequency coordinates are derived as `(0:nfrequency-1) * frequency_spacing`.
Sample count is derived from storage, so reconstruction with different row counts
changes duration/extent. Constructors validate intrinsic metadata, one-based
storage axes, and channel count. No coordinates or channel values are type parameters.

## Array operations and selection

Same-shape arithmetic with ordinary arrays preserves metadata. Multiple series
operands must agree on coordinates, epoch, and ordered channels.
Broadcasting validates the operands and destination before executing on parent
storage. Use `flatten` around repeated numerical work to avoid repeated checks.
Shape-changing arithmetic, irregular sample selections, reductions, reshape,
transpose, and dimension permutations return ordinary arrays.

```julia
2 .* x                         # TimeSeries, unchanged metadata
x[:, :H1]                      # ordinary vector
x[:, [:H1]]                    # one-channel TimeSeries
x[:, [:L1, :H1]]                # channels in the requested order
x[101:200, :]                   # copy, advanced epoch
view(x, 101:200, :)             # TimeSeries wrapping shared storage
view(x, 1:2:8191, :)            # ordinary view
```

Unknown channel names throw. Repeated or empty channel selections and dimension-dropping
selections return ordinary arrays/views. Only contiguous increasing time slices preserve
TimeSeries metadata. Strided, reversed, or irregular selections return ordinary arrays.
Frequency selections preserve FrequencySeries only when all rows remain in order;
partial bands return ordinary arrays. Channel-only selection preserves either series type.
`crop(x; start, stop)` and `pad(x; before, after)` apply only to TimeSeries.
Crop uses grid-aligned half-open boundaries; padding before shifts the epoch.

`frequency_bins(h; fmin, fmax)` returns indices of existing bins inside inclusive bounds,
clipped to the available grid, without rounding bounds onto a new grid. Disjoint bounds
produce an empty range; nonfinite or reversed bounds throw. Keep these indices separately
when selecting a band.

Coordinate comparison retains exact fast paths and a few-ulp tolerance: four ulps
of absolute coordinates plus eight ulps of displacement and spacing. Ambiguous
comparisons whose tolerance reaches half a sample/bin are rejected. BigFloat
metadata is available for GPS epochs requiring sub-microsecond precision.

## Physical AbstractFFTs transforms

GWSeries uses the AbstractFFTs API and does not load an FFT backend. Load one
supporting your parent storage (for example FFTW for CPU arrays):

```julia
using AbstractFFTs, FFTW
h = rfft(x)                    # FrequencySeries: dt*rfft(parent(x), 1)
z = irfft(h, nsamples(x))       # TimeSeries, original epoch and channels
h = to_frequencydomain(x)      # real samples only
z = to_timedomain(h; nsamples=8192)
```

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
inverse = plan_to_timedomain(h; nsamples=nsamples(x))
metadata(forward)               # output constructor keywords
transform!(h, forward, x)
raw_h = transform(forward, parent(x))
```

Prepared transforms use AbstractFFTs plans, physical scaling, and reusable private
scratch storage. Inputs are preserved. Plans are not reentrant; use separate
plans for concurrent tasks. Raw execution accepts prevalidated arrays; wrapped
execution also checks metadata. Backend-specific preparation keywords pass to
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
this contract. StaticArrays interoperability uses a small optional StaticArraysCore extension.
ArrayInterface integration is optional and does not advertise
pointer/stride capabilities or change host-versus-traced kernel selection.

## Inverting a selected band

A band is an ordinary array, not an independently invertible spectrum. Restore its bins
explicitly to produce a band-limited signal; zeros do not recover discarded values:

```julia
bins = frequency_bins(h; fmin=20.1, fmax=1024.0)
band = view(h, bins, :)
full = similar(h)
fill!(full, 0)
full[bins, :] = band
t_bandlimited = irfft(full, nsamples(x))
```

No original time length is stored in FrequencySeries. Changing the number of bins while
holding Δf fixed changes the time sampling implied by a chosen inverse length.
