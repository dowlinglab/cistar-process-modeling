import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.regenerate_archived_figures import (  # noqa: E402
    OUTPUT_NAMES,
    c4_plus_percent,
    region_checkpoint,
    rok_checkpoint,
)


class ArchivedFigureRegenerationTests(unittest.TestCase):
    def test_exactly_eight_publication_artifacts_are_covered(self):
        self.assertEqual(set(OUTPUT_NAMES), {"F3", "F5", "F6", "F7", "F8", "FS1", "FS2", "FS5"})

    def test_checkpoint_names_preserve_m5_sequential_suffix(self):
        self.assertTrue(rok_checkpoint("M5").name.endswith("_sequential_solve.json.gz"))
        self.assertFalse("sequential" in rok_checkpoint("M2").name)
        self.assertTrue(region_checkpoint("EF-Basin").is_file())

    def test_c4_plus_percentage_uses_butene_through_nonene(self):
        flows = {"pentane": 10.0, "butene": 2.0, "pentene": 3.0, "nonene": 5.0}
        self.assertAlmostEqual(c4_plus_percent(flows), 50.0)


if __name__ == "__main__":
    unittest.main()
