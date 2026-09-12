import sys
import unittest
from pathlib import Path

from scipy import sparse


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.fingerprint_m5_bakken_derivatives import (  # noqa: E402
    _parse_perturbation,
    canonical_sparse_fingerprints,
)


class DerivativeFingerprintTests(unittest.TestCase):
    def test_sparse_fingerprint_is_storage_order_independent(self):
        first = sparse.coo_matrix(
            ([2.0, -3.0], ([0, 1], [1, 0])), shape=(2, 2)
        )
        second = sparse.coo_matrix(
            ([-3.0, 2.0], ([1, 0], [0, 1])), shape=(2, 2)
        )
        names = ["a", "b"]
        self.assertEqual(
            canonical_sparse_fingerprints(first, names, names),
            canonical_sparse_fingerprints(second, names, names),
        )

    def test_perturbation_parser(self):
        self.assertEqual(
            _parse_perturbation("H105_temperature=1.25"),
            ("H105_temperature", 1.25),
        )


if __name__ == "__main__":
    unittest.main()
