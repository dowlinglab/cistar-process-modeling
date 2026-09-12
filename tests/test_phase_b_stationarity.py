import sys
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.diagnose_m5_bakken_stationarity import (  # noqa: E402
    reduced_space_gradient,
    stationarity_projection,
)


class StationarityProjectionTests(unittest.TestCase):
    def test_gradient_in_row_span_has_negligible_residual(self):
        rows = sparse.csr_matrix([[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        result = stationarity_projection(np.array([2.0, 2.0, -1.0]), rows)
        self.assertLess(result["stationarity_residual_norm_2"], 1e-9)

    def test_projection_retains_null_space_component(self):
        rows = sparse.csr_matrix([[1.0, 1.0]])
        result = stationarity_projection(np.array([1.0, -1.0]), rows)
        self.assertAlmostEqual(
            result["stationarity_residual_norm_2"], np.sqrt(2.0), places=9
        )
        self.assertAlmostEqual(result["relative_stationarity_residual"], 1.0)

    def test_reduced_gradient_eliminates_state_equation(self):
        # c = state + 2*design = 0 and f = 3*state + design, so
        # df/ddesign = 1 + 3*(-2) = -5.
        result = reduced_space_gradient(
            np.array([3.0, 1.0]),
            sparse.csr_matrix([[1.0, 2.0]]),
            [1],
        )
        self.assertAlmostEqual(result["reduced_objective_gradient"][0], -5.0)


if __name__ == "__main__":
    unittest.main()
