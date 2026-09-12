import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from regenerate_composite_curves import (  # noqa: E402
    CurveCase,
    add_isothermal_duties,
    cooling_utility_mw,
    curve_coordinates,
    discover_cases,
    parse_case,
)


class CompositeCurveRegenerationTests(unittest.TestCase):
    def test_discovers_every_archived_optimal_curve(self):
        cases = discover_cases(REPO_ROOT / "results" / "solution_data.xlsx")
        self.assertEqual(len(cases), 23)
        self.assertEqual(len({case.filename for case in cases}), 23)

    def test_normalizes_misspelled_basin_sheet(self):
        case = parse_case("HI_M5_tax=0.045_EF-Basn_optimal")
        self.assertIsNotNone(case)
        self.assertEqual(case.region, "EF-Basin")
        self.assertIn("EF-Basin", case.filename)

    def test_tax_case_uses_reported_tonne_units(self):
        case = CurveCase(
            sheet="test", model_code=5, tax_token="0.045",
            tax_usd_per_kg=0.045, region="Bakken"
        )
        tables = {
            "taxes": pd.DataFrame({"C-tax-rate": [45.0], "Qw": [22.528]}),
            "models": pd.DataFrame(),
            "regions": pd.DataFrame(),
        }
        self.assertEqual(cooling_utility_mw(case, tables), 22.528)

    def test_isothermal_duty_keeps_coordinates_finite(self):
        frame = pd.DataFrame(
            [
                [0, "T inlet", "K", 300.0, 600.0, 700.0],
                [1, "T outlet", "K", 500.0, 600.0, 350.0],
                [2, "Heat duty", "W", 1.0e6, 2.0e6, -3.0e6],
            ]
        )
        coordinates = curve_coordinates(frame, cooling_utility=1.0)
        for values in coordinates.values():
            self.assertTrue(np.isfinite(values).all())
        cold_t = np.asarray(coordinates["cold_temperature_k"])
        cold_q = np.asarray(coordinates["cold_mw"])
        at_reactor = np.flatnonzero(cold_t == 600.0)
        self.assertEqual(len(at_reactor), 2)
        self.assertAlmostEqual(cold_q[at_reactor[1]] - cold_q[at_reactor[0]], 2.0)

    def test_isothermal_jump_shifts_only_higher_temperatures(self):
        temperature, heat = add_isothermal_duties(
            np.array([300.0, 500.0, 700.0]),
            np.array([0.0, 10.0, 20.0]),
            [(500.0, 3.0)],
        )
        np.testing.assert_allclose(temperature, [300.0, 500.0, 500.0, 700.0])
        np.testing.assert_allclose(heat, [0.0, 10.0, 13.0, 23.0])


if __name__ == "__main__":
    unittest.main()
