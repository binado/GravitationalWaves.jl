"""Select CI matrix entries and verify the aggregate result (standard library only)."""

import json
import os
import sys

PACKAGES = (
    "AbstractGWWaveforms",
    "LALWaveforms",
    "GWParameterConversion",
    "GWWaveforms",
    "GWDistributions",
    "GWDetectors",
    "GWTimeFrequencyDomains",
    "GWNoise",
    "GWMatchedFiltering",
    "GWSurrogates",
)
AFFECTED = {
    "GWTimeFrequencyDomains": {"GWTimeFrequencyDomains", "GWNoise", "GWDetectors", "GWMatchedFiltering"},
    "AbstractGWWaveforms": {"AbstractGWWaveforms", "LALWaveforms", "GWWaveforms", "GWSurrogates"},
    "LALWaveforms": {"LALWaveforms", "GWWaveforms"},
    "GWParameterConversion": {"GWParameterConversion", "LALWaveforms", "GWWaveforms"},
    "GWWaveforms": {"GWWaveforms"},
    "GWDistributions": {"GWDistributions"},
    "GWDetectors": {"GWDetectors"},
    "GWNoise": {"GWNoise"},
    "GWMatchedFiltering": {"GWMatchedFiltering"},
    "GWSurrogates": {"GWSurrogates"},
}
ROOT_DOCUMENTATION = {"README.md", "AGENTS.md", "CLAUDE.md", "LICENSE"}


def select_packages(paths):
    if not isinstance(paths, list) or any(
        not isinstance(path, str) or not path for path in paths
    ):
        raise ValueError("Changed files must be a JSON array of nonempty paths")
    selected = set()
    for path in paths:
        if path in ROOT_DOCUMENTATION:
            continue
        parts = path.split("/")
        if len(parts) >= 3 and parts[0] == "packages" and parts[1] in AFFECTED:
            selected.update(AFFECTED[parts[1]])
        elif path.startswith("python/surrogates/") or path == "python/generate_surrogates.py":
            selected.update(AFFECTED["GWSurrogates"])
        elif path.startswith("python/"):
            selected.update(AFFECTED["GWParameterConversion"] | {"GWDetectors", "GWMatchedFiltering"})
        else:
            # Shared configuration, new packages and unfamiliar inputs run everything.
            selected.update(PACKAGES)
    return [package for package in PACKAGES if package in selected]


def select_matrix(event, changed_files, changed_count=0):
    if event == "pull_request":
        packages = select_packages(json.loads(changed_files))
        # GitHub's PR files API returns at most 3,000 files. Never use a partial diff.
        if changed_count > 3000:
            packages = list(PACKAGES)
        return {"packages": packages, "julia": ["1.12"]}
    if event in {"push", "workflow_dispatch"}:
        return {"packages": list(PACKAGES), "julia": ["1.12", "1"]}
    raise ValueError(f"Unsupported CI event: {event}")


def check_status(needs):
    for job in ("changes", "format"):
        if needs[job]["result"] != "success":
            raise ValueError(f"{job} did not succeed: {needs[job]['result']}")
    packages = json.loads(needs["changes"]["outputs"]["packages"])
    if (
        not isinstance(packages, list)
        or any(not isinstance(p, str) or p not in PACKAGES for p in packages)
        or len(packages) != len(set(packages))
    ):
        raise ValueError("Invalid selected package array")
    expected = "success" if packages else "skipped"
    if needs["test"]["result"] != expected:
        raise ValueError(f"Expected tests to be {expected}: {needs['test']['result']}")


def main():
    if sys.argv[1:] == ["select"]:
        matrix = select_matrix(
            os.environ["CI_EVENT"],
            os.environ.get("CHANGED_FILES", ""),
            int(os.environ.get("PR_CHANGED_FILE_COUNT", "0")),
        )
        with open(os.environ["GITHUB_OUTPUT"], "a") as output:
            output.writelines(
                f"{name}={json.dumps(values)}\n" for name, values in matrix.items()
            )
        summary = (
            f"Packages: {', '.join(matrix['packages']) or 'none (documentation only)'}\n\n"
            f"Julia: {', '.join(matrix['julia'])}\n"
        )
        print(summary)
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as output:
            output.write(summary)
    elif sys.argv[1:] == ["check"]:
        check_status(json.loads(os.environ["CI_NEEDS"]))
        print("All required CI jobs passed.")
    else:
        raise ValueError("Usage: ci.py select|check")


if __name__ == "__main__":
    main()
