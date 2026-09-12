import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class PhaseACoverageTests(unittest.TestCase):
    def test_every_manifest_item_has_one_disposition(self):
        manifest = json.loads(
            (REPO_ROOT / "reproducibility" / "result_manifest.json").read_text()
        )
        coverage = json.loads(
            (REPO_ROOT / "reproducibility" / "phase_a_coverage.json").read_text()
        )
        expected = [
            item["id"]
            for section in ("main_figures", "supporting_figures", "supporting_tables")
            for item in manifest[section]
        ]
        observed = [item["id"] for item in coverage["items"]]
        self.assertEqual(len(observed), 22)
        self.assertEqual(observed, expected)
        self.assertEqual(len(observed), len(set(observed)))
        self.assertTrue(coverage["pr_gate"]["ready"])

    def test_gate_does_not_hide_partial_or_nonreproduced_cases(self):
        coverage = json.loads(
            (REPO_ROOT / "reproducibility" / "phase_a_coverage.json").read_text()
        )
        combined = " ".join(
            item["status"] + " " + item["assessment"] for item in coverage["items"]
        )
        for required in ("EF-2", "EF-8", "EF-10", "EF-11", "IPOPT"):
            self.assertIn(required, combined + " " + coverage["historical_environment"]["limitation"])


if __name__ == "__main__":
    unittest.main()
