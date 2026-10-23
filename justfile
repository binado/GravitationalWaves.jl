# Override with JULIA='julia +1.12' to choose a Julia version.
julia := env_var_or_default("JULIA", "julia")

default:
    @just --list

# Resolve the workspace (root + all members) into a single manifest.
resolve:
    {{ julia }} --startup-file=no --project=. -e 'using Pkg; Pkg.resolve()'

# Run Pkg.test() for every workspace member.
test:
    for pkg in packages/AbstractWaveforms packages/LALWaveforms packages/GWParameterConversion packages/FrequencyDomainWaveforms packages/GWDistributions packages/GWDetectors packages/GWSeries packages/GWNoise packages/GWMatchedFiltering packages/SurrogateWaveforms packages/AbstractWaveformsTests; do {{ julia }} --startup-file=no --project=$pkg -e 'using Pkg; Pkg.test()'; done

# Run Pkg.test() for a single workspace member, e.g. `just test-pkg LALWaveforms`.
test-pkg pkg:
    {{ julia }} --startup-file=no --project=packages/{{ pkg }} -e 'using Pkg; Pkg.test()'

# Install the exact versions recorded in the shared workspace manifest.
instantiate:
    {{ julia }} --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'

# Bump versions across the whole workspace into a new shared manifest.
update:
    {{ julia }} --startup-file=no --project=. -e 'using Pkg; Pkg.update()'

# Format everything in place (CI enforces Runic).
fmt:
    runic --inplace .

# Generate bilby reference fixtures (gitignored) and print the LAL libdir export.
fixtures:
    uv sync --project python
    uv run --project python python python/generate.py
    uv run --project python python python/generate_detectors.py
    uv run --project python python python/generate_inference.py
    uv sync --project python/surrogates
    uv run --project python/surrogates python python/generate_surrogates.py --data-dir reference/surrogates-data --download
    @echo 'export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)'
    @echo 'export GWSURROGATES_DATA_DIR='$(pwd)/reference/surrogates-data

# Run all heavy AD integration suites locally (no LAL or Python setup needed).
test-integration: test-integration-series test-integration-waveforms test-integration-detectors test-integration-inference test-integration-surrogates test-reference-waveforms test-reference-surrogates

# FrequencyDomainWaveforms against LALWaveforms as the reference generator; local only, needs LALWAVEFORMS_LIBDIR.
test-reference-waveforms:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/FrequencyDomainWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/FrequencyDomainWaveforms/test_reference.jl")'

# ForwardDiff and Mooncake waveform gradients against finite differences.
test-integration-waveforms:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/FrequencyDomainWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/FrequencyDomainWaveforms/runtests.jl")'

# Geometry, timing, and phase-sensitive projection gradients (ForwardDiff and Mooncake).
test-integration-detectors:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWDetectors -e 'using Pkg; Pkg.instantiate(); include("integration/GWDetectors/runtests.jl")'

# Joint metric gradients and waveform/projection/likelihood composition.
test-integration-inference:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWMatchedFiltering -e 'using Pkg; Pkg.instantiate(); include("integration/GWMatchedFiltering/runtests.jl")'

# NRSur7dq4/v2 AD: ForwardDiff and Mooncake. Downloads the model data on first use (or set GWSURROGATES_DATA_DIR); the compile times are several minutes.
test-integration-surrogates:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/SurrogateWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/SurrogateWaveforms/runtests.jl")'

# Layout comparison; output includes timings, allocations, CPU, and Julia version.
benchmark-inference:
    {{ julia }} --startup-file=no --project=benchmarks/GWMatchedFiltering -e 'using Pkg; Pkg.instantiate()'
    {{ julia }} --startup-file=no --threads=1 --project=benchmarks/GWMatchedFiltering benchmarks/GWMatchedFiltering/layouts.jl

# Physical Fourier AD against explicit DFTs (ForwardDiff and Mooncake).
test-integration-series:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWSeries -e 'using Pkg; Pkg.instantiate(); include("integration/GWSeries/runtests.jl")'

# Native surrogates against upstream gwsurrogate (PyCall reference generator); local only.
test-reference-surrogates:
    uv sync --project python/surrogates --locked
    PYTHON=$(uv run --project python/surrogates python -c 'import sys; print(sys.executable)') {{ julia }} --project=integration/SurrogateWaveforms/gwsurrogate_reference -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
    {{ julia }} --check-bounds=yes --project=integration/SurrogateWaveforms/gwsurrogate_reference integration/SurrogateWaveforms/gwsurrogate_reference/runtests.jl
