# GWTimeFrequencyDomains

Detached immutable metadata for uniformly sampled data. Domains do not contain
samples or implement array interfaces. Numerical operations accept ordinary
arrays, with vectors for one channel and matrices in `(sample, channel)` layout.

```julia
using GWTimeFrequencyDomains

td = TimeDomain(; epoch=1126259460.125, nsamples=8192,
                sampling_frequency=4096, channels=(:H1,))
fd = frequency_domain(td)
x = zeros(Float32, nsamples(td))
check_data(td, x)                     # once, on the host
```

`FrequencyDomain(; frequencies=0:0.5:2048, epoch=t0, channels=(:H1,))`
constructs a frequency domain directly, including offset bands. Epochs are
explicit seconds (GPS for detector data); channels are ordered unique symbols,
with strings canonicalized to symbols. Constructors validate intrinsic metadata.
Metadata precision is at least Float64, preserving BigFloat inputs independently
of the arrays' element types. Epoch and spacing values are fields, not type parameters.

## Accessors and explicit checks

`epoch`, `channels`, and `nsamples` apply to both domains. The latter counts
frequency bins for frequency domains. `times(td)` and `frequencies(fd)` are lazy
ranges. `sample_spacing`, `sampling_frequency`, and `duration = N*dt` apply only to
time domains; `frequency_spacing` applies to frequency domains.

`check_data(domain, array)` checks numeric vector/matrix storage, one-based axes,
row count, and channel count. A one-channel domain accepts vectors and one-column
matrices. `check_same_grid(a,b)` checks domains, grid counts/spacings/origins,
epochs, and channels. `check_same_frequency_grid(fa,fb)` omits epoch/channel checks.
Consumers implement `check_compatible(domain, prepared_object)` for stationary
noise/metrics, whose compatibility deliberately ignores epoch.

Grid comparisons accept identical coordinates exactly, otherwise use four ulps
of each absolute coordinate plus eight ulps of relative displacement and spacing.
Ambiguous differences whose tolerance reaches half a sample/bin are rejected.
Subtracting epochs before scaling preserves GPS-scale timing precision. Frequency
epoch comparisons use the inferred even FFT spacing, not an observed original
sampling frequency. Keep each operand's domain until checking alignment; detached
arrays cannot establish where they came from.

`crop(domain, values; start, stop)` copies a half-open interval with on-grid bounds.
`pad(domain, values; before=0, after=0)` zero-pads counts. Both validate on the host
and return `(new_domain, new_values)`. Padding before shifts time epoch or frequency
origin; negative frequencies are rejected. Empty and singleton domains keep their
origins and spacings. These operations are outside the compiled inference contract.

## Physical Fourier conversion

```julia
using FFTW
check_fourier_domain(td)
h = to_frequencydomain(td, x)  # dt*rfft(x, 1); arrays only
x2 = to_timedomain(fd, h)      # irfft(h, N, 1)/dt

forward = plan_to_frequencydomain(td, x)
inverse = plan_to_timedomain(fd, h)
output_domain(forward)
h2 = transform(forward, x)
transform!(h2, forward, x)     # separate output storage, input unchanged
```

`frequency_domain(td)` and `time_domain(fd)` derive and validate metadata once.
Forward conversion requires even `N >= 2`; inverse conversion requires a zero-origin
grid with `M >= 2` and infers `N=2(M-1)`, `dt=1/(N*df)`. The phase origin is relative
to the epoch; there is no absolute-GPS rephasing. Imaginary DC/Nyquist parts are
discarded. Removing low bins prevents inversion until restored explicitly; removing
high bins lowers inverse sampling frequency while preserving `1/df` duration.

Plans validate prototypes on the host and cache backend plans, scaling, shape/type
requirements and output domains. Execution assumes matching prevalidated arrays;
backend checks may remain. `transform!` requires separate, nonaliasing destination
storage of the plan's output shape/type. Inverse plans own reusable scratch input
to preserve source values and support real-valued spectra. Plans are not reentrant;
use separate plans per concurrent task. CPU Float32/Float64 arrays and strided views
are supported, with matching complex frequency storage and unchanged sample precision.
There is no automatic host conversion for unsupported backend arrays.

Domain-first numerical methods perform no repeated metadata validation. Fixed domains
can be captured in compiled inference functions, while physical parameters and samples
remain runtime inputs. Dynamic whole-domain Reactant arguments and FFT differentiation
are not established by this API.
