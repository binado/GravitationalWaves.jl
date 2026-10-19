# LALSuite binary release preparation

This directory contains the candidate Yggdrasil recipe for `LALSuite_jll`.
**It is not a released binary provider.** LALWaveforms still requires its existing
external-library configuration. Adoption is blocked until the JLL is registered
and all four required platforms pass binary audits and native validation.

Upstream draft: [Yggdrasil #15187](https://github.com/JuliaPackaging/Yggdrasil/pull/15187).

## Source and build contract

- Upstream: LALSuite 7.26.16, matching `python/uv.lock`; LAL 7.7.1 and
  LALSimulation 6.2.1.
- Protected tag: [`lalsuite-v7.26.16`](https://git.ligo.org/lscsoft/lalsuite/-/tags/lalsuite-v7.26.16),
  commit `548bc5710bea021239da56e9beff1feaae163956`.
- Archive SHA-256: `ca263715e3378577b90103ce4c7de78acccac4afecb87ba44b85807e497964b1`
  (downloaded and checked locally; complete archive is 240,334,392 bytes).
- Targets: `x86_64-linux-gnu`, `aarch64-linux-gnu`, `x86_64-apple-darwin`,
  `aarch64-apple-darwin`. Windows and musl are deferred.
- Shared LAL, LALSupport, and LALSimulation, headers, pkg-config metadata, and
  library licenses. No supported executable products, SWIG bindings, Python
  embedding, documentation, OpenMP, fast-math, or host-specific optimizations.
- Build dependency versions: FFTW_jll 3.3.11, GSL_jll 2.8.1, HDF5_jll 1.12.2,
  Zlib_jll 1.3.1; GCC 8 toolchain. The serial HDF5 release avoids requiring MPI
  and has artifacts for all four targets. Dependency build revisions and all
  transitive dependencies must be recorded from actual builds before release.

LALSimulation probes Python embedding independently of `--disable-python`.
The recipe explicitly disables that probe before generating configure scripts.
Every LAL library product uses `dont_dlopen=true`, so LALWaveforms can select an
external provider before loading any bundled dependencies.

HDF5 is enabled for users supplying ROM/surrogate datasets through `LAL_DATA_PATH`.
The recipe does not download or bundle those datasets. No fallback data directory
is compiled into the libraries.

## Build and validate

Copy `build_tarballs.jl` and `bundled/` into `L/LALSuite/` in a Yggdrasil checkout. Use a Linux
BinaryBuilder environment, or a supported Docker runner on macOS. First build
Linux x86-64 and macOS ARM64, then the other two targets:

```sh
cd L/LALSuite
julia --project=/path/to/binarybuilder build_tarballs.jl x86_64-linux-gnu,aarch64-apple-darwin --verbose
julia --project=/path/to/binarybuilder build_tarballs.jl aarch64-linux-gnu,x86_64-apple-darwin --verbose
```

Keep BinaryBuilder's relocation/dependency auditor enabled. Do not use
`--skip-audit`, `--ignore-audit-errors`, or `--skip-dlopen` as release evidence.
The `dont_dlopen` products require additional native loading checks.

Generate a candidate wrapper using BinaryBuilder's local JLL deployment flow.
On each native target, create a temporary Julia 1.12 environment containing that
wrapper and the local AbstractGWWaveforms, GWParameterConversion, and LALWaveforms
packages (`Pkg.develop` them together). Relocate the artifacts outside their build
tree before testing. Run:

```sh
env -u LALWAVEFORMS_LIBDIR julia --startup-file=no --project=/path/to/validation-env packaging/LALSuite/smoke.jl
```

`smoke.jl` verifies import/JLL laziness, every required symbol, cached first use,
TaylorF2, IMRPhenomD, precessing IMRPhenomXPHM, transformed TaylorT4, and absence
of loaded libpython. It uses the existing explicit provider to test candidate
artifacts before adoption; it does not test the future default-provider branch.

For each native build, retain `ldd -r` (Linux) or `otool -L` and `otool -l` (macOS)
output for **all three libraries and their transitive dependencies**. Require no
unresolved dependencies, Python linkage, references to build directories, or
absolute paths to the build machine. Check installed pkg-config files too, and
confirm headers/licenses and HDF5 support. Record compressed/unpacked artifact
sizes, dependency versions/build revisions, compiler flags, and audit logs.

Generate bilby fixtures independently using the locked `python/` environment.
Run complete LALWaveforms and GWWaveforms suites against the candidate artifacts
in separate Julia processes. Keep existing tolerances; record failures and
investigate compiler/dependency/libm differences before considering changes.

## Release gates and follow-up adoption

| Target | Build / relocation audit | Native smoke | Both complete suites | Artifact sizes |
| --- | --- | --- | --- | --- |
| Linux x86-64 glibc | Pending | Pending | Pending | Pending |
| macOS ARM64 | Pending | Pending | Pending | Pending |
| Linux ARM64 glibc | Pending | Pending | Pending | Pending |
| macOS x86-64 | Pending | Pending | Pending | Pending |

Local preparation verified the source archive, BinaryBuilder recipe metadata on
Julia 1.12 and 1.13, patch application, recipe shell syntax, Runic formatting, and
the 14 CI helper tests. The existing external loader tests and all four waveform
smoke inputs passed against wheel libraries on Julia 1.12. These checks do not
validate the candidate artifacts. The macOS ARM64 build attempt stopped before
compilation because Docker is unavailable; all native artifact results are pending.

Submit the recipe upstream, complete all four validation rows, and wait for
General registration before adopting the JLL. An existing upstream
[core-LAL proposal](https://github.com/JuliaPackaging/Yggdrasil/pull/14202)
does not supply LALSimulation; coordinate packaging with its maintainers.

The subsequent repository PR must:

1. Add the registered JLL UUID and restrict compatibility to the validated release.
2. Resolve nonempty environment override, then saved preference, then the JLL.
   Load the JLL with `Base.require` only inside the bundled branch, using explicit
   product paths and the existing lock/global-LAL-before-LALSimulation flow.
   Cache state only after symbol resolution and error-handler configuration.
3. Preserve external discovery, immediate setter loading, restart warnings,
   explicit-override failures, and actionable unsupported-platform errors.
4. Test all provider choices in isolated subprocesses, including precompilation,
   import laziness, invalid overrides, precedence, and restart behavior.
5. Make normal LALWaveforms/GWWaveforms CI use the JLL; retain bilby fixture
   generation and a separate Linux wheel-provider job. Add four Julia 1.12 native
   installation jobs with no Python setup, preferences, or overrides; include
   them in required status aggregation and its helper tests.
6. Update ordinary installation/setup documentation and development recipes only
   once the default provider exists. Retain Python fixture and external-provider
   instructions, and explain `LAL_DATA_PATH` without automatic dataset downloads.

Do not merge a temporary `[sources]` path, private artifact URL, fabricated UUID,
or a reduced target matrix. Future upstream upgrades repeat all validation gates
and remain separate from packaging changes. Numerical/AD kernels are untouched;
heavy AD integration suites are unnecessary for this preparation.
