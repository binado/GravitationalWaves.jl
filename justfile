# Override with JULIA='julia +1.12' to choose a Julia version.
julia := env_var_or_default("JULIA", "julia")

default:
    @just --list

# Resolve the workspace (root + all members) into a single manifest.
resolve:
    {{ julia }} --startup-file=no --project=. -e 'using Pkg; Pkg.resolve()'

# Run Pkg.test() for every workspace member.
test:
    for pkg in packages/AbstractGWWaveforms packages/LALWaveforms packages/GWParameterConversion packages/GWWaveforms packages/GWDistributions packages/GWDetectors packages/GWNoise packages/GWMatchedFiltering; do {{ julia }} --startup-file=no --project=$pkg -e 'using Pkg; Pkg.test()'; done

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
    @echo 'export LALWAVEFORMS_LIBDIR=$(uv run --project python python python/libdir.py)'

# Run all heavy AD integration suites locally (no LAL or Python setup needed).
test-integration: test-integration-conversion test-integration-waveforms test-integration-detectors test-integration-inference

# Reactant CPU compilation and Enzyme gradients for parameter conversions.
test-integration-conversion:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWParameterConversion -e 'using Pkg; Pkg.instantiate(); include("integration/GWParameterConversion/runtests.jl")'

# ForwardDiff, Enzyme and Mooncake waveform gradients against finite differences, and Reactant CPU TaylorF2.
test-integration-waveforms:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/GWWaveforms/runtests.jl")'

# Native and compiled geometry, timing, and phase-sensitive projection gradients.
test-integration-detectors:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWDetectors -e 'using Pkg; Pkg.instantiate(); include("integration/GWDetectors/runtests.jl")'

# Joint metric gradients and waveform/projection/likelihood composition.
test-integration-inference:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWMatchedFiltering -e 'using Pkg; Pkg.instantiate(); include("integration/GWMatchedFiltering/runtests.jl")'

# Layout comparison; output includes timings, allocations, CPU, and Julia version.
benchmark-inference:
    {{ julia }} --startup-file=no --project=benchmarks/GWMatchedFiltering -e 'using Pkg; Pkg.instantiate()'
    {{ julia }} --startup-file=no --threads=1 --project=benchmarks/GWMatchedFiltering benchmarks/GWMatchedFiltering/layouts.jl
