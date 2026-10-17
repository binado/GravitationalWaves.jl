# Override with JULIA='julia +1.12' to choose a Julia version.
julia := env_var_or_default("JULIA", "julia")

default:
    @just --list

# Resolve the workspace (root + all members) into a single manifest.
resolve:
    {{ julia }} --startup-file=no --project=. -e 'using Pkg; Pkg.resolve()'

# Run both heavy AD integration suites locally (no LAL or Python setup needed).
test-integration: test-integration-conversion test-integration-waveforms

# Reactant CPU compilation and Enzyme gradients for parameter conversions.
test-integration-conversion:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWParameterConversion -e 'using Pkg; Pkg.instantiate(); include("integration/GWParameterConversion/runtests.jl")'

# ForwardDiff, Enzyme and Mooncake waveform gradients against finite differences.
test-integration-waveforms:
    {{ julia }} --startup-file=no --check-bounds=yes --project=integration/GWWaveforms -e 'using Pkg; Pkg.instantiate(); include("integration/GWWaveforms/runtests.jl")'
