import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.compare_m5_scaling_policies import (  # noqa: E402
    compare_factor_maps,
)


class ScalingPolicyComparisonTests(unittest.TestCase):
    def test_factor_ratios_and_missing_values(self):
        result = compare_factor_maps(
            {"same": 2.0, "ratio": 1e-3, "missing": None},
            {"same": 2.0, "ratio": 1e-1, "missing": 1.0},
            top=3,
        )
        self.assertEqual(result["reference_missing_factors"], 1)
        self.assertEqual(result["candidate_missing_factors"], 0)
        self.assertEqual(result["counts_above_absolute_log10_ratio"]["1.0"], 2)
        self.assertEqual(result["largest_factor_ratios"][0]["name"], "missing")


if __name__ == "__main__":
    unittest.main()
