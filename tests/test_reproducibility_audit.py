import json
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from audit_published_tables import audit_snapshot  # noqa: E402


class PublishedReferenceTests(unittest.TestCase):
    def test_reference_and_manifest_are_valid_json(self):
        for relative in (
            "reproducibility/published_reference.json",
            "reproducibility/result_manifest.json",
        ):
            with (REPO_ROOT / relative).open(encoding="utf-8") as stream:
                self.assertEqual(json.load(stream)["schema_version"], 1)

        run_records = sorted((REPO_ROOT / "reproducibility" / "runs").glob("*.json"))
        self.assertGreaterEqual(len(run_records), 2)
        for path in run_records:
            with path.open(encoding="utf-8") as stream:
                record = json.load(stream)
            self.assertEqual(record["schema_version"], 1, path.name)
            self.assertIn("experiment_id", record, path.name)

    def test_manifest_covers_every_published_figure_and_table(self):
        with (REPO_ROOT / "reproducibility/result_manifest.json").open(
            encoding="utf-8"
        ) as stream:
            manifest = json.load(stream)
        self.assertEqual(
            [item["id"] for item in manifest["main_figures"]],
            [f"F{i}" for i in range(1, 9)],
        )
        self.assertEqual(
            [item["id"] for item in manifest["supporting_figures"]],
            [f"FS{i}" for i in range(1, 9)],
        )
        self.assertEqual(
            [item["id"] for item in manifest["supporting_tables"]],
            [f"TS{i}" for i in range(1, 7)],
        )

    def test_postprocessed_snapshot_has_only_classified_mismatches(self):
        report = audit_snapshot("postprocessed")
        self.assertEqual(report["comparison_count"], 350)
        self.assertEqual(report["missing_row_count"], 0)
        self.assertEqual(report["unexpected_mismatch_count"], 0)
        self.assertEqual(report["known_mismatch_count"], 47)
        issue_counts = {}
        for item in report["mismatches"]:
            issue = item["known_issue"]
            issue_counts[issue] = issue_counts.get(issue, 0) + 1
        self.assertEqual(issue_counts["EMISSIONS-NORMALIZATION-001"], 25)
        self.assertEqual(issue_counts["PUBLISHED-SNAPSHOT-DRIFT-001"], 22)

    def test_migrated_snapshot_records_pre_postprocessing_economics(self):
        report = audit_snapshot("migrated")
        self.assertEqual(report["comparison_count"], 350)
        self.assertEqual(report["missing_row_count"], 0)
        self.assertEqual(report["unexpected_mismatch_count"], 0)
        issue_counts = {}
        for item in report["mismatches"]:
            issue = item["known_issue"]
            issue_counts[issue] = issue_counts.get(issue, 0) + 1
        self.assertEqual(
            issue_counts["LEGACY-ECONOMICS-POSTPROCESSING-001"], 29
        )


if __name__ == "__main__":
    unittest.main()
