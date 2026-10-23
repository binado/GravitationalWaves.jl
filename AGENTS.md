# Repository agent guidance

This repository is a Julia monorepo for gravitational-wave modeling.
The root is a workspace, not an umbrella package.

## Layout and scope

- Independent packages live under `packages/`, each with its own API, version,
  dependencies, tests, and license.
- Root `Project.toml` defines workspace membership; do not assume every package
  directory belongs to the workspace or the default CI matrix.
- Read a package's local `AGENTS.md` before changing it, including when changing
  several packages together.
- Use the root [README](README.md) and package READMEs for setup and detailed
  behavior. Keep agent guidance focused on constraints for future edits.

## julia styleguide

- Avoid the `@inbounds` macro
- Avoid closures
- Prefer map / broadcasting / comprehensions over index-matched inputs when you can
- Prefer dedicated indexing methods like `eachindex`, `pairs`, `axes`,
`firstindex`, `lastindex` over ranges
- Use `BenchmarkTools` for benchmarking performance
- When designing a public API, prefer function names and signatures that mirror
julia stdlib/Base/ecossystem conventions (dispatch on new types but avoid type pyracy)

## Commands

Run these from the repository root; see the [justfile](justfile) for all recipes
and their prerequisites. Set `JULIA` to select a Julia executable or version.

```sh
just --list                 # available recipes
just instantiate            # install recorded workspace dependencies
just resolve                # resolve the shared workspace manifest
just update                 # update workspace dependencies
just test-pkg <Package>     # test one package
just test                   # test all workspace members
just fmt                    # format with Runic
```

Use dependency-changing recipes only when needed by the task. Preserve unrelated
manifest and dependency changes already present in the working tree.

## Verification

- Use `Pkg.test()` (directly or through the recipes) so test-only dependencies
  are resolved. Directly including a test file may miss its test environment.
- Test affected packages and downstream consumers when changing shared behavior.
  Consult [.github/ci/ci.py](.github/ci/ci.py) for package-selection rules rather
  than maintaining a second dependency map here.
- Run relevant integration suites for numerical kernels or AD changes before
  submitting them. These use standalone environments under `integration/`;
  ordinary package tests do not run them.
- Check [.github/workflows/CI.yml](.github/workflows/CI.yml) for current Julia
  versions, setup, formatting, and required checks. Local integration acceptance
  remains necessary even when default CI passes.
- Run `git diff --check` and appropriate formatting before finishing. Pure
  documentation edits need no Julia tests unless they change executable examples.

## Git and review

- Use Conventional Commits, work on a non-`main` branch, and open PRs for changes.
- Keep diffs focused on the task and preserve unrelated work.
- Describe the resulting behavior, validation performed, and material limitations
  in the PR; link detailed documentation instead of copying it into guidance.
