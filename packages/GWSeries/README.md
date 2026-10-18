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
sample spacing; this does not give the container an observed duration. There is
no implicit alignment, rephasing, resampling, or channel reordering.

Host construction validates metadata. Sample element and storage types remain
generic. Validate/unwrap on the host before calling compiled or differentiated
array kernels and reconstruct results with `with_data` outside compilation.
Whole-series Reactant arguments and FFT differentiation are not supported contracts.
