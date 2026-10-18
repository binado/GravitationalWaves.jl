# GWSeries

Explicit time/frequency sampling containers. Series do not subtype `AbstractArray`
and do not provide indexing, arithmetic, or broadcasting. Numerical work uses
`parent(series)` and existing array kernels.

```julia
using GWSeries
x = TimeSeries(zeros(Float32, 4096); epoch=1126259460.125,
               sampling_frequency=4096, channels=(:H1,))
y = with_data(x, ones(Float32, 4096))
check_same_grid(x, y)
segment = crop(x; start=epoch(x), stop=epoch(x) + 0.5)
padded = pad(segment; before=128, after=128)
```

Vectors require one channel label; matrices use `(sample, channel)` layout and
require one unique label per column. Strings canonicalize to symbols. Epochs are
explicit timestamps in seconds (GPS for detector data). Metadata uses at least
Float64 precision, preserving BigFloat inputs, regardless of sample precision.
`parent` returns the original array. Mutating it is visible; resizing it is
unsupported. `copy` copies samples; `with_data` retains replacement storage and
requires identical complete shape. One-based axes are required.

`nsamples` counts rows, including frequency bins. `times` and `frequencies` return
lazy ranges. `sample_spacing`, `sampling_frequency`, and half-open segment
`duration = nsamples * sample_spacing` apply only to time series.

`crop` copies the half-open interval `[start, stop)`; both boundaries must align
with the grid and lie within its half-open extent (the endpoint is accepted as a
boundary). Empty and singleton series retain their origin and spacing. `pad`
zero-fills sample/bin counts; padding before shifts the time epoch or frequency
origin. Frequency padding cannot create negative frequencies.

Alignment subtracts origins before scaling. Its tolerance is four ulps of each
absolute coordinate plus eight ulps of the relative displacement and spacing.
Coordinates whose tolerance reaches half a sample/bin are ambiguous and rejected.
`check_same_grid` checks counts, spacing, accumulated drift, origins, epochs, and
ordered channels. Frequency-series epoch comparison uses the inferred even FFT
sample spacing; this does not give the container an observed duration. Cropping and padding are host operations; compiled support is not claimed. There is
no implicit alignment, rephasing, resampling, or channel reordering.

Host construction validates metadata. Sample element and storage types remain
generic. Validate/unwrap on the host before calling compiled or differentiated
array kernels and reconstruct results with `with_data` outside compilation.
Whole-series Reactant arguments and FFT differentiation are not supported contracts.

## Fourier conversion

Load `FFTW` to activate CPU Float32/Float64 array and array-view conversion:

```julia
using FFTW
h = to_frequencyseries(x) # dt * rfft(parent(x), 1), even N >= 2
x2 = to_timeseries(h)     # irfft(parent(h), N, 1) / dt
```

Transforms act along the first dimension and preserve epoch/channel labels. The
phase origin is relative to the epoch; no absolute-GPS phase factor is applied.
Inverse conversion requires at least two bins starting at exactly zero and
infers `N = 2(M-1)`, `dt = 1/(N*df)`. Imaginary DC/Nyquist components are discarded
by `irfft`. Odd time series may be stored but must be explicitly cropped/padded
before forward conversion. Removing low bins prevents inversion until explicitly
restored (for example with `pad`). Removing high bins lowers the inverse sampling
frequency while keeping duration `1/df`; pad back to restore the earlier rate.
Conversion allocates output storage and retains FFTW sample precision (Float32 or
Float64); metadata precision remains independent. Unsupported devices/AD element types do not
fall back to host materialization.
