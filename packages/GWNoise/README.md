# GWNoise

Read one-sided PSD/ASD tables and draw stationary Gaussian strain noise in Julia.
PSD estimation and continuous streaming are outside the current scope.

```julia
using GWNoise, GWTimeFrequencyDomains, Random

curve = read_psd("H1_psd.txt")  # or read_asd("H1_asd.txt")
S = psd(range(20, 1024; step = 0.25), curve)

td = TimeDomain(;epoch=1126259460.0, nsamples=8192, sampling_frequency=2048, channels=(:H1,))
model = GaussianNoise(
    curve, td;
    minimum_frequency = 20.0, maximum_frequency = 1024.0,
)
rng = Xoshiro(42)
f = frequencies(model)
n_fd = noise_fd(rng, model)
n_td = noise_td(rng, model)  # a new, independent draw
```

Files have two columns: frequency in Hz and PSD in strain²/Hz (or ASD in strain/√Hz).
Whitespace or comma delimiters, blank lines and `#` comments are supported. Readers
accept a path or `IO`. There is no inferred header, unit conversion or network access.
`TabulatedPSD(f, values; kind=:psd)` also accepts in-memory tables. Tables and draws
use `Float64`/`ComplexF64`.

The default interpolation is linear in PSD and frequency. Use
`read_psd(path; interpolation=:loglog)` for linear interpolation of log PSD against
log frequency; this requires positive frequencies and densities. ASD is squared
before interpolation. `asd(f, curve)` evaluates the square root of the resulting PSD.
Evaluation outside file coverage throws; `psd(f, curve; outside=:zero)` explicitly
requests zero fill. Frequencies must be strictly increasing; densities nonnegative.

`GaussianNoise` prepares the uniform FFT grid and effective PSD once. Its default
band is the file coverage intersected with `[0, sampling_frequency/2]`; explicit
bands must be covered by the file. Draw scales are computed per bin when sampling.
Other bins, DC, and even-length Nyquist are zero.
Odd-length final bins are sampled normally. `psd(model)` returns the effective PSD
including these zeros. No automatic taper or smoothing changes the target.

For `N = nsamples` and `fs = sampling_frequency`, duration is `N/fs`, frequency
spacing is `fs/N`, and the FD vector has `N÷2+1` elements. Interior real and imaginary
parts independently have variance `PSD/(4Δf)`. The physical Fourier convention is:

```julia
n_td = to_timedomain(model.time_domain, n_fd)
```

Each draw is a finite periodic FFT segment. Concatenating independent segments does
not produce a continuous stream. Use the transform above for matching FD and TD
realizations; successive sampler calls draw independent noise. Omit the RNG to use
`Random.default_rng()`. RNGs are never reseeded by this package.

`noise_fd!([rng], out, model)` fills an existing complex vector. `noise_td`
allocates its output and FFT temporaries; neither FFT plans nor mutable scratch
buffers are cached in the model.

Run tests with `julia --project=packages/GWNoise -e 'using Pkg; Pkg.test()'` from
the monorepo root. No Python or LAL installation is needed.

## Domain preparation and explicit validation

`GaussianNoise(curve, td; ...)` takes a single-channel `TimeDomain` and prepares
stationary numerical state. Epoch does not affect a draw. Draws return arrays;
`noise_fd` and `noise_td` consume only their prepared model. There are no sample
wrappers or series-returning overloads.

Before repeated draws into `out`, call `check_data(domain, out)`
and `check_compatible(domain, model)` explicitly. Compatibility checks sampling
and channel identity, ignoring epoch; numerical draw methods do not repeat these
checks. Invalid output shapes are outside their preconditions. `frequencies`
is shared with GWTimeFrequencyDomains.

All FFT conversions are implemented by GWTimeFrequencyDomains; GWNoise loads FFTW
and LinearAlgebra to activate that package's FFT extension. Odd time-domain lengths
remain supported through `to_timedomain(model.time_domain, n_fd)`, which uses the
known target length. `frequency_domain(td)` still requires even lengths.

Version 0.2 replaces `GaussianNoise(curve; sampling_frequency, nsamples, ...)`
with `GaussianNoise(curve, td; ...)` and requires separate output validation.
