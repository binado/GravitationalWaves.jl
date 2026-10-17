# AbstractGWWaveforms

The interface shared by the gravitational-wave waveform backends in this repository,
[`GWWaveforms`](../GWWaveforms) (pure Julia, differentiable) and [`LALWaveforms`](../LALWaveforms)
(LALSimulation through `ccall`). It has no dependencies.

It defines `AbstractWaveformGenerator` and the generic function

```julia
hphc_fd(f, params, generator) -> (; plus, cross)
```

where `f` is a uniform frequency grid starting at 0, `params` a `NamedTuple` with bilby's
parameter names, and `generator` a concrete generator such as `TaylorF2()`, `IMRPhenomD()` or
`LALWaveformGenerator("IMRPhenomXPHM")`. Code written against this package works with any backend.

and `validate_parameters(f, params, generator)`, which throws `ArgumentError` on invalid input and
returns `nothing` otherwise. Backends that compilers can trace (`GWWaveforms` with Reactant) do not
check parameter values in `hphc_fd`; callers run `validate_parameters` once on the host instead.
The fallback accepts everything, so generic code can always call it.

Every generator provides `approximant`, `reference_frequency`, `minimum_frequency` and
`maximum_frequency`, and may answer the capability queries `spin_support` (`:none`, `:aligned` or
`:precessing`), `supports_mode_array` and `is_differentiable`. `check_frequency_grid` and
`frequency_bounds` hold the grid checks and masking bounds backends share.

`test/contract.jl` defines `test_generator_contract(g; params, f)`, which backends `include` from
their tests to check a generator against the contract.

Licensed under GPL-2.0-or-later.
