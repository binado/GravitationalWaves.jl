# SurrogateWaveforms integration tests

The native-vs-gwsurrogate comparison lives in [`gwsurrogate_reference/`](gwsurrogate_reference/runtests.jl): the
PyCall adapter [`SurrogatePyCallWaveforms`](../../packages/SurrogatePyCallWaveforms/README.md)
is the reference generator for `AbstractWaveformsTests.test_reference`. It covers both models'
signed modes (model grid, fixed time step, explicit times) and windowed polarizations by
mismatch. Run it locally with `just test-reference-surrogates`; it needs the pinned Python
environment (`uv sync --project python/surrogates --locked`) and model data from the shared
cache or `GWSURROGATES_DATA_DIR`. CI does not run it.

`just test-integration-surrogates` runs the ForwardDiff and Mooncake integration tests and needs no Python.
