# GWNoise

Read one-sided PSD/ASD tables and draw stationary Gaussian strain noise in Julia.
PSD estimation and continuous streaming are outside the current scope.

```julia
using GWNoise, GWSeries, Random, FFTW

curve = read_psd("H1_psd.txt")  # or read_asd("H1_asd.txt")
S = psd.(range(20, 1024; step = 0.25), curve)

td = TimeSeries(zeros(8192); epoch=1126259460.0, sampling_frequency=2048, channels=(:H1,))
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
Use `psd.(frequencies, curve)` or `asd.(frequencies, curve)` for collections of
frequencies; each curve broadcasts as a scalar.
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
n_td = to_timedomain(n_fd; nsamples=nsamples(td))
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

## Series and metadata preparation

`GaussianNoise(curve, prototype::TimeSeries; ...)` prepares a single-channel
stationary segment. Alternatively supply `nsamples`, `epoch`, `channels`, and
`sampling_frequency` or `sampling_interval` as keywords without allocating a prototype.
Models store sampling metadata and the effective PSD, independently of sample buffers.

Draws return one-channel FrequencySeries or TimeSeries. `noise_fd!` accepts either
series or ordinary vector/matrix storage; series compatibility checks grid and channels
while ignoring epoch. Raw-buffer methods check shape. The original time sample count
is retained for correct odd/even inverse transforms.

GWSeries implements physical transforms through AbstractFFTs. Load a backend such as
FFTW when drawing time-domain noise; GWNoise does not load a backend implicitly.
Frequency-domain draws and PSD preparation do not need an FFT backend.

## Bundled detector noise curves

`available_psds()` lists twelve bundled designs: A+, Virgo O5 low/high, KAGRA
128 Mpc, ET D and COBA 10/15/20 km, and CE reference, 20 km,
20 km post-merger, and 40 km low-frequency curves. `psd_path(name)` accepts
a symbol or string and returns a local table path.

```julia
using GWNoise
curve = read_psd(psd_path(:AplusDesign))
value = psd(100.0, curve)
```

The supplied files are retained verbatim in `data/noise_curves`; columns are
frequency in Hz and one-sided PSD in strain²/Hz. No additional provenance or
licensing information was supplied with the tables. Interpolation and coverage
follow `read_psd` and `psd`, including errors outside the tabulated band by default.

GWDetectors provides independent geometry and default sensitivity mappings:

```julia
using GWDetectors
curve = read_psd(psd_path(sensitivities.H1))
```
