import sys
import unittest
from pathlib import Path

from pyomo.environ import ConcreteModel, Constraint, Var


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproducibility"))

from diagnose_m5_region_state import (  # noqa: E402
    _badly_scaled_variables,
    _largest_residuals,
    _near_bound_variables,
)


class RegionStateDiagnosticTests(unittest.TestCase):
    def test_largest_residuals_are_ranked(self):
        model = ConcreteModel()
        model.x = Var(initialize=1.0)
        model.small = Constraint(expr=model.x == 1.5)
        model.large = Constraint(expr=model.x == 3.0)

        records = _largest_residuals(model, tolerance=0.1, limit=1)

        self.assertEqual(records, [{"name": "large", "absolute_residual": 2.0}])

    def test_bad_scaling_and_near_bounds_are_serializable(self):
        model = ConcreteModel()
        model.large = Var(initialize=1e6)
        model.small = Var(initialize=1e-6)
        model.bounded = Var(bounds=(0.0, 10.0), initialize=0.0)

        scaling_names = {
            record["name"] for record in _badly_scaled_variables(model, limit=10)
        }
        near_bounds = _near_bound_variables(model, limit=10)

        self.assertEqual(scaling_names, {"large", "small"})
        self.assertEqual(
            near_bounds,
            [
                {
                    "name": "bounded",
                    "value": 0.0,
                    "lower_bound": 0.0,
                    "upper_bound": 10.0,
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
