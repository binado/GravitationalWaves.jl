# Format with Runic (https://github.com/fredrikekre/Runic.jl), installed in the shared
# `@runic` environment so it never touches the package's own dependencies.

default:
    @just --list

# One-time: install Runic into the shared @runic environment.
install-runic:
    julia --startup-file=no --project=@runic -e 'using Pkg; Pkg.add("Runic")'

# Format every Julia file in place.
fmt:
    julia --startup-file=no --project=@runic -m Runic --inplace .

# Fail if any file is not formatted (what CI runs).
fmt-check:
    julia --startup-file=no --project=@runic -m Runic --check .

# Run the test suite on Julia 1.12 (the supported version).
test:
    julia +1.12 --project -e 'using Pkg; Pkg.test()'
