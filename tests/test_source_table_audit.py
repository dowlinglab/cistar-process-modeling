import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from audit_source_tables import audit_source_tables  # noqa: E402


class SourceTableAuditTests(unittest.TestCase):
    def setUp(self):
        self.report = audit_source_tables()

    def test_all_source_table_rows_are_classified(self):
        self.assertEqual(self.report["finding_count"], 33)
        self.assertEqual(
            self.report["status_counts"],
            {
                "match": 28,
                "explained_source_split": 2,
                "missing_public_source": 2,
                "implementation_discrepancy": 1,
            },
        )
        self.assertEqual(self.report["unexpected_count"], 0)

    def test_h104_bounds_are_not_silently_treated_as_matching(self):
        finding = next(
            item
            for item in self.report["findings"]
            if item["table"] == "table_s1" and item["name"] == "T_H104"
        )
        self.assertEqual(finding["status"], "implementation_discrepancy")
        self.assertEqual(finding["implemented"], {"lower": 273.15, "upper": 1500.0})

    def test_public_source_omissions_remain_explicit(self):
        missing = {
            item["name"]
            for item in self.report["findings"]
            if item["status"] == "missing_public_source"
        }
        self.assertEqual(missing, {"methane_recovery_fraction", "methane_gwp"})


if __name__ == "__main__":
    unittest.main()
