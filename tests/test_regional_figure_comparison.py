import sys
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from compare_regional_figure_data import (  # noqa: E402
    _compare_frames,
    _gen_curves,
    _workbook_region,
)


class RegionalFigureComparisonTests(unittest.TestCase):
    def test_numeric_and_unit_differences_are_classified_separately(self):
        fresh = pd.DataFrame(
            [["mol/s", 2.5], ["K", "-"]],
            index=["flow", "temperature"],
            columns=["Units", "s01"],
        )
        archived = pd.DataFrame(
            [["mole / second", 2.0], ["K", "-"]],
            index=["flow", "temperature"],
            columns=["Units", "s01"],
        )

        result = _compare_frames(fresh, archived)

        self.assertEqual(result["numeric_cells_compared"], 1)
        self.assertEqual(result["maximum_absolute_difference"], 0.5)
        self.assertEqual(result["text_mismatch_count"], 0)
        self.assertEqual(result["unit_label_mismatch_count"], 1)

    def test_ef_basin_heat_sheet_preserves_archived_typo(self):
        self.assertEqual(_workbook_region("EF-Basin"), "EF-Basn")
        self.assertEqual(_workbook_region("EF-7"), "EF-7")

    def test_curve_generation_matches_nearest_neighbor_interpolation(self):
        temperature, heat = _gen_curves(
            pd.Series([300.0]).to_numpy(),
            pd.Series([400.0]).to_numpy(),
            pd.Series([50.0]).to_numpy(),
        )
        self.assertEqual(temperature.tolist(), [300.0, 400.0])
        self.assertEqual(heat.tolist(), [50.0, 0.0])


if __name__ == "__main__":
    unittest.main()
