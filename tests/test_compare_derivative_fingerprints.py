import sys
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.compare_derivative_fingerprints import (  # noqa: E402
    _compare,
    _matrix_entries,
    _vector_entries,
)


class CompareDerivativeFingerprintTests(unittest.TestCase):
    def test_matrix_entries_align_names(self):
        first = sparse.coo_matrix(([2.0], ([0], [1])), shape=(2, 2))
        second = sparse.coo_matrix(([2.0], ([1], [0])), shape=(2, 2))
        first_entries = _matrix_entries(first, np.array(["a", "b"]), np.array(["a", "b"]))
        second_entries = _matrix_entries(second, np.array(["b", "a"]), np.array(["b", "a"]))
        self.assertEqual(first_entries, second_entries)

    def test_symmetric_entries_ignore_triangle_orientation(self):
        upper = sparse.coo_matrix(([3.0], ([0], [1])), shape=(2, 2))
        lower = sparse.coo_matrix(([3.0], ([1], [0])), shape=(2, 2))
        names = np.array(["a", "b"])
        self.assertEqual(
            _matrix_entries(upper, names, names, symmetric=True),
            _matrix_entries(lower, names, names, symmetric=True),
        )

    def test_compare_reports_scaled_difference(self):
        result = _compare({("a",): 100.0}, {("a",): 100.01}, top=1)
        self.assertAlmostEqual(result["maximum_scaled_difference"], 0.01 / 100.01)

    def test_vector_entries_omit_exact_zeros(self):
        self.assertEqual(
            _vector_entries(np.array([0.0, 2.0]), np.array(["a", "b"])),
            {("b",): 2.0},
        )


if __name__ == "__main__":
    unittest.main()
