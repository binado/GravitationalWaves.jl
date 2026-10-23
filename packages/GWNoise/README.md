# GWNoise

Read one-sided PSD/ASD tables and draw stationary Gaussian strain noise in Julia.
PSD estimation and continuous streaming are outside the current scope.

```julia
using GWNoise, Random

curve = read_psd("H1_psd.txt")  # or read_asd("H1_asd.txt")
f = range(0.0; step = 0.25, length = 4097)          # one-sided grid from DC
S = reshape(band_psd(curve, f; fmin = 20.0, fmax = 1024.0), :, 1)  # (nf, nc)
W = whitening_factor(S)                              # for GWMatchedFiltering
n = sample_noise(Xoshiro(42), S, step(f))            # (nf, nc) ComplexF64 draw
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

## Plain-array conventions

Everything works on ordinary arrays; no series types are involved.

- Data are `(nf, nc)` complex one-sided arrays from DC, in `rfft(x) * dt` units.
- `S` is the one-sided PSD in strain²/Hz: `(nf, nc)` real for independent channels, or
  `(nc, nc, nf)` Hermitian for correlated channels. Zeros mark inactive bins; arrays always
  have the full length `nf` and bands are zeros, never slices.
- Coefficients of a draw satisfy `E[n n†] = S / (2Δf)`, so interior real and imaginary
  parts independently have variance `S/(4Δf)`.
- `band_psd(curve, f; fmin, fmax)` evaluates a `TabulatedPSD` on a frequency range and
  zeroes everything outside `[fmin, fmax]`, plus DC and the last bin (it assumes an
  even-length transform, so for an odd-length one it also drops a real bin).
- `whitening_factor(S)` returns `W = L⁻¹` with `L L† = S`: `1/√S` for `(nf, nc)` input
  (zero where `S` is zero) or the inverse lower Cholesky factor per bin for `(nc, nc, nf)`
  input. Blocks are validated once (finite, Hermitian, positive definite) and invalid
  ones throw `ArgumentError`.
- `sample_noise([rng], S, df)` and `sample_noise!([rng], out, S, df)` draw
  coefficients with that covariance. Random numbers are drawn for every bin, active or not,
  so a stream stays aligned across bands. DC and Nyquist are not checked; they should be
  zero in `S` by convention.

Each draw is a finite periodic FFT segment; concatenating independent segments does not
produce a continuous stream. RNGs are never reseeded by this package. Transforming to the
time domain is left to the caller (for instance with GWSeries and an FFT backend).

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
