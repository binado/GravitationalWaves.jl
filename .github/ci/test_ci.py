import json
import unittest

from ci import PACKAGES, check_status, select_matrix, select_packages


CONVERSION_AFFECTED = ["LALWaveforms", "GWParameterConversion", "FrequencyDomainWaveforms"]


class SelectionTests(unittest.TestCase):
    def test_package_dependencies(self):
        cases = {
            "GWSeries": ["AbstractWaveforms", "LALWaveforms", "FrequencyDomainWaveforms", "GWDetectors", "GWSeries", "GWNoise", "GWMatchedFiltering", "SurrogateWaveforms"],
            "GWDistributions": ["GWDistributions"],
            "GWDetectors": ["GWDetectors"],
            "GWNoise": ["GWNoise"],
            "GWMatchedFiltering": ["GWMatchedFiltering"],
            "SurrogateWaveforms": ["SurrogateWaveforms"],
            "FrequencyDomainWaveforms": ["FrequencyDomainWaveforms"],
            "LALWaveforms": ["LALWaveforms", "FrequencyDomainWaveforms"],
            "AbstractWaveforms": ["AbstractWaveforms", "LALWaveforms", "FrequencyDomainWaveforms", "SurrogateWaveforms"],
            "GWParameterConversion": CONVERSION_AFFECTED,
        }
        for package, expected in cases.items():
            for suffix in (
                "src/model.jl",
                "test/deleted.jl",
                "Project.toml",
                "README.md",
            ):
                with self.subTest(package=package, suffix=suffix):
                    self.assertEqual(
                        select_packages([f"packages/{package}/{suffix}"]), expected
                    )

    def test_shared_and_unknown_paths(self):
        for path in (
            "Project.toml",
            "Manifest.toml",
            "Manifest-v1.12.toml",
            ".github/workflows/CI.yml",
            ".github/ci/ci.py",
            ".gitignore",
            "packages/NewPackage/src/new.jl",
            "new-tool.sh",
            "docs/new.md",
        ):
            with self.subTest(path=path):
                self.assertEqual(select_packages([path]), list(PACKAGES))

    def test_python_tooling(self):
        for name in (
            "generate.py",
            "generate_detectors.py",
            "generate_inference.py",
            "libdir.py",
            "cases.json",
            "uv.lock",
            "pyproject.toml",
        ):
            with self.subTest(name=name):
                self.assertEqual(
                    select_packages([f"python/{name}"]),
                    CONVERSION_AFFECTED + ["GWDetectors", "GWMatchedFiltering"],
                )

    def test_surrogate_python_tooling_selects_only_surrogates(self):
        for path in (
            "python/generate_surrogates.py",
            "python/surrogates/uv.lock",
            "python/surrogates/pyproject.toml",
        ):
            with self.subTest(path=path):
                self.assertEqual(select_packages([path]), ["SurrogateWaveforms"])

    def test_documentation_only_and_empty(self):
        self.assertEqual(
            select_packages(["README.md", "AGENTS.md", "CLAUDE.md", "LICENSE"]), []
        )
        self.assertEqual(select_packages([]), [])

    def test_union_and_deduplication(self):
        self.assertEqual(
            select_packages(
                [
                    "packages/GWDistributions/test/new.jl",
                    "packages/LALWaveforms/src/fd.jl",
                    "packages/LALWaveforms/Project.toml",
                    "README.md",
                ]
            ),
            ["LALWaveforms", "FrequencyDomainWaveforms", "GWDistributions"],
        )

    def test_shared_waveform_fixtures(self):
        self.assertEqual(select_packages(["test/waveforms/parameters.csv"]),
                         ["FrequencyDomainWaveforms", "SurrogateWaveforms"])
        for path in ("test/waveforms/surrogate_oracle/runtests.jl",
                     "integration/SurrogatePyCallWaveforms/Project.toml"):
            self.assertEqual(select_packages([path]), ["SurrogateWaveforms"])

    def test_cross_package_rename(self):
        # paths-filter supplies both the deleted old path and added new path.
        self.assertEqual(
            select_packages(
                [
                    "packages/GWDistributions/src/old.jl",
                    "packages/GWParameterConversion/src/new.jl",
                ]
            ),
            CONVERSION_AFFECTED + ["GWDistributions"],
        )

    def test_unusual_filenames_are_data(self):
        self.assertEqual(
            select_packages(["packages/FrequencyDomainWaveforms/test/space ' quote $(cmd)\n.jl"]),
            ["FrequencyDomainWaveforms"],
        )

    def test_event_matrices(self):
        self.assertEqual(
            select_matrix("pull_request", '["README.md"]'),
            {"packages": [], "julia": ["1.12"]},
        )
        for event in ("push", "workflow_dispatch"):
            self.assertEqual(
                select_matrix(event, ""),
                {"packages": list(PACKAGES), "julia": ["1.12", "1"]},
            )

    def test_large_pr_falls_back_to_full_matrix(self):
        self.assertEqual(
            select_matrix("pull_request", '["README.md"]', 3001),
            {"packages": list(PACKAGES), "julia": ["1.12"]},
        )

    def test_invalid_input_fails(self):
        for data in ("", "null", "{}", '"README.md"', "[null]", '[""]', "[1]"):
            with self.subTest(data=data), self.assertRaises(ValueError):
                select_matrix("pull_request", data)
        with self.assertRaises(ValueError):
            select_matrix("unknown", "[]")


class StatusTests(unittest.TestCase):
    @staticmethod
    def needs(packages, test="success", changes="success", formatting="success"):
        return {
            "changes": {
                "result": changes,
                "outputs": {"packages": json.dumps(packages)},
            },
            "test": {"result": test},
            "format": {"result": formatting},
            "surrogate-oracle": {"result": "success" if isinstance(packages, list) and "SurrogateWaveforms" in packages else "skipped"},
        }

    def test_oracle_failure_blocks_status(self):
        for result in ("failure", "cancelled", "skipped"):
            needs = self.needs(["SurrogateWaveforms"])
            needs["surrogate-oracle"]["result"] = result
            with self.assertRaises(ValueError):
                check_status(needs)

    def test_success(self):
        check_status(self.needs(list(PACKAGES)))
        check_status(self.needs(["FrequencyDomainWaveforms"]))
        check_status(self.needs([], test="skipped"))

    def test_unsuccessful_dependencies(self):
        for result in ("failure", "cancelled", "skipped"):
            for field in ("changes", "formatting", "test"):
                with (
                    self.subTest(result=result, field=field),
                    self.assertRaises(ValueError),
                ):
                    check_status(self.needs(list(PACKAGES), **{field: result}))

    def test_empty_selection_requires_skipped_tests(self):
        for result in ("success", "failure", "cancelled"):
            with self.subTest(result=result), self.assertRaises(ValueError):
                check_status(self.needs([], test=result))

    def test_invalid_outputs_fail(self):
        for packages in (
            None,
            {},
            "",
            ["unknown"],
            [1],
            [["FrequencyDomainWaveforms"]],
            ["FrequencyDomainWaveforms", "FrequencyDomainWaveforms"],
        ):
            with self.subTest(packages=packages), self.assertRaises(ValueError):
                check_status(self.needs(packages))
        needs = self.needs([])
        needs["changes"]["outputs"] = {}
        with self.assertRaises(KeyError):
            check_status(needs)


if __name__ == "__main__":
    unittest.main()
