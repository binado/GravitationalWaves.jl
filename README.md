# LALWaveforms.jl

Frequency-domain compact-binary waveforms from [LALSimulation](https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/),
called directly through `ccall`, with [bilby](https://bilby-dev.github.io/bilby/)'s parameter names, units and
parameter-to-LAL conversions.

```julia
using LALWaveforms

settings = WaveformSettings("IMRPhenomXPHM"; reference_frequency = 50.0, minimum_frequency = 20.0)
f = range(0, 1024; step = 0.25)
params = (
    mass_1 = 36.0, mass_2 = 29.0,          # detector-frame solar masses
    luminosity_distance = 500.0,           # Mpc
    theta_jn = 0.4, phase = 1.3,
    a_1 = 0.5, a_2 = 0.3, tilt_1 = 1.0, tilt_2 = 0.6, phi_12 = 0.3, phi_jl = 4.0,
)
(; plus, cross) = hphc_fd(f, params, settings)
```

`WaveformSettings` builds the LAL waveform dictionary once; reuse it across calls. Optional parameters
(`a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, eccentricity`) default to 0, and unknown keys throw.

Scope for now: frequency domain only, uniform `AbstractRange` grids starting at 0, detector-frame masses,
no tidal parameters.

## Pointing the package at LAL

LAL is not bundled. Provide the directory holding `liblal` and `liblalsimulation` in one of two ways
(the environment variable wins):

- `LALWAVEFORMS_LIBDIR=/path/to/libs`
- `LALWaveforms.set_libdir!("/path/to/libs")`, saved as a preference, or
  `LALWaveforms.set_libdir_from_python!()` to take it from the `lalsuite` wheel in the active Python environment.

The `lalsuite` wheels carry the libraries in `lal/.dylibs/` (macOS) or `lalsuite.libs/` (Linux).
Loading is lazy: the package loads without LAL and the first waveform call explains how to configure it.

## Tests

The reference tests compare against bilby output generated locally (not committed):

```sh
uv sync --project test/python
uv run --project test/python python test/python/generate.py
export LALWAVEFORMS_LIBDIR=$(uv run --project test/python python test/python/libdir.py)
julia --project -e 'using Pkg; Pkg.test()'
```

Without fixtures the reference tests warn and skip; with `CI=true` they fail.

## Formatting

Code is formatted with [Runic](https://github.com/fredrikekre/Runic.jl):

```sh
julia -e 'using Pkg; Pkg.Apps.add("Runic")'
runic --inplace .
```

## License

GPL-2.0-or-later.
