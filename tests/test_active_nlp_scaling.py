import sys
import unittest
from pathlib import Path

from pyomo.environ import Block, ConcreteModel, Suffix, Var


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.run_m5_bakken_tax_series import (  # noqa: E402
    DESIGN_VARIABLE_NAMES,
    FileDeterminism,
    _clear_scaling_suffixes,
)


class ActiveNlpScalingTests(unittest.TestCase):
    def test_published_design_variable_names_are_stable(self):
        self.assertEqual(len(DESIGN_VARIABLE_NAMES), 8)
        self.assertEqual(len(set(DESIGN_VARIABLE_NAMES)), 8)

    def test_file_determinism_symbol_sort_is_available(self):
        self.assertEqual(FileDeterminism.SORT_SYMBOLS.name, "SORT_SYMBOLS")

    def test_clear_scaling_suffixes_recurses_and_preserves_other_suffixes(self):
        model = ConcreteModel()
        model.x = Var(initialize=2.0)
        model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.scaling_factor[model.x] = 0.5
        model.other = Suffix(direction=Suffix.EXPORT)
        model.other[model.x] = 3.0
        model.block = Block()
        model.block.y = Var(initialize=4.0)
        model.block.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.block.scaling_factor[model.block.y] = 0.25

        summary = _clear_scaling_suffixes(model)

        self.assertEqual(summary, {"suffixes_cleared": 2, "entries_cleared": 2})
        self.assertEqual(len(model.scaling_factor), 0)
        self.assertEqual(len(model.block.scaling_factor), 0)
        self.assertEqual(model.other[model.x], 3.0)


if __name__ == "__main__":
    unittest.main()
