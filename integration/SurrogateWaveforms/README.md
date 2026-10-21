# SurrogateWaveforms integration tests

`test_pycall.jl` compares the native Julia generator with the upstream Python oracle
through `SurrogatePyCallWaveforms`. It covers both models' signed modes and polarizations
on the model grid, a fixed time step, and explicit sample times, including a nontrivial
initial frame. Model data uses the shared cache or `GWSURROGATES_DATA_DIR`.

Configure Python following [the adapter setup](../SurrogatePyCallWaveforms/README.md), then run
the oracle comparisons from the repository root:

```sh
julia --project=integration/SurrogateWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/SurrogateWaveforms/test_pycall.jl")'
```

`just test-integration-surrogates` runs these comparisons along with the ForwardDiff
and Mooncake integration tests. The Python setup is required for that full suite too.
