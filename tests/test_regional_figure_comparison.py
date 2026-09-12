import sys
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from compare_regional_figure_data import (  # noqa: E402
    _compare_frames,
    _figure_7_classification,
    _gen_curves,
    _workbook_region,
)
from compare_bakken_figure_data import _case_sheet, _run_tax_rate  # noqa: E402
from regenerate_composite_curves import CurveCase  # noqa: E402
from src.result_extraction import dataframe_payload as _dataframe_payload  # noqa: E402


class RegionalFigureComparisonTests(unittest.TestCase):
    def test_bakken_sheet_lookup_uses_recorded_tax_token(self):
        cases = [CurveCase("HI_M5_tax=1e-05_Bakken_optimal", 5, "1e-05", 1e-5, "Bakken")]
        self.assertEqual(_case_sheet(cases, 5, 1e-5).tax_token, "1e-05")
        with self.assertRaises(ValueError):
            _case_sheet(cases, 3, 1e-5)

    def test_bakken_tax_rate_works_for_series_and_single_case_records(self):
        self.assertEqual(_run_tax_rate({"case": {}}, {"co2_tax_usd_per_kg": 0.045}), 0.045)
        self.assertEqual(_run_tax_rate({"case": {"co2_tax_usd_per_kg": 0.045}}, {}), 0.045)

    def test_dataframe_payload_preserves_trace_float_precision(self):
        trace_value = 1.260558903831854e-14
        payload = _dataframe_payload(
            pd.DataFrame([[trace_value]], index=["trace"], columns=["s01"])
        )
        self.assertEqual(payload["data"][0][0], trace_value)

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
        self.assertEqual(result["maximum_material_relative_difference"], 0.25)
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

    def test_figure_7_snapshot_drift_requires_workbook_agreement(self):
        self.assertEqual(
            _figure_7_classification(-1.2, 3e-7),
            "PUBLISHED-SNAPSHOT-DRIFT-001",
        )
        self.assertEqual(
            _figure_7_classification(-1.2, 2e-2),
            "unclassified_fresh_solution_difference",
        )


if __name__ == "__main__":
    unittest.main()
