# SurrogatePyCallWaveforms

Experimental adapter that calls upstream `gwsurrogate` through PyCall and implements the
`AbstractWaveforms` time-domain mode contract, so it can serve as a reference generator for
`AbstractWaveformsTests.test_reference`. MIT licensed, like `gwsurrogate`. It is kept out of the
root workspace and CI package matrix, since PyCall needs a configured Python: its tests and the
native-vs-oracle comparison run locally only.

## Build PyCall with uv's interpreter

The root uv project keeps the optional reference dependencies in the `gwsurrogate` extra.
From the repository root, create or update that environment:

```sh
uv sync --extra gwsurrogate
```

PyCall chooses and records its Python interpreter when it is built. Set `PYTHON` to the
executable inside the uv environment **before** instantiating/building PyCall. The commands below
work from the repository root and can be repeated after changing the uv environment or Julia
depot:

```sh
export PYTHON="$(uv run --extra gwsurrogate python -c 'import sys; print(sys.executable)')"
uv run --extra gwsurrogate python -c 'import gwsurrogate, numpy; print(gwsurrogate.__file__)'
julia --project=packages/SurrogatePyCallWaveforms -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
julia --project=packages/SurrogatePyCallWaveforms -e 'using PyCall; println(PyCall.python)'
```

The last command should print the same interpreter path shown by `uv run ... sys.executable`.
Setting `PYTHON` only for a later Julia run will not switch an existing PyCall build; run
`Pkg.build("PyCall")` again with the desired path. PyCall loads that interpreter and its installed
packages for the lifetime of the Julia process.

## Load and use a model

Use the existing HDF5 file from the Julia cache or provide a local data directory; this adapter
does not download a second copy:

```julia
using SurrogateWaveforms, SurrogatePyCallWaveforms
path = SurrogateWaveforms.data_path(:NRSur7dq4)
g = SurrogatePyCallWaveforms.load_generator(path; model=:NRSur7dq4, ell_max=4, f_ref=25.0)
```

The resulting generator accepts `total_mass`, `luminosity_distance`, `mass_ratio`, `spin_1`,
`spin_2`, `inclination`, and `phi_ref` parameters, with optional `init_orbphase` and `init_quat`.
`waveform_modes_td` returns the contract's ordered signed labels and `(sample, mode)` matrix in
physical strain. `hphc_td` is inherited from `AbstractWaveforms` and uses the LAL sky-angle
convention `(inclination, π/2 - phi_ref)`. Python evaluation is not differentiable.

The adapter's own tests (`julia --project=packages/SurrogatePyCallWaveforms -e 'using Pkg; Pkg.test()'`,
after the PyCall build above) check the contract, units and argument validation. The native
implementation is compared with it by `just test-reference-surrogates`, which runs
`integration/SurrogateWaveforms/gwsurrogate_reference/runtests.jl`.
