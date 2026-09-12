import io
import sys
import unittest
from pathlib import Path

import idaes


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

IS_IDAES_212 = idaes.__version__.startswith("2.12")
if IS_IDAES_212:
    from reproducibility.diagnose_phase_b_m5 import _capture_report  # noqa: E402


@unittest.skipUnless(IS_IDAES_212, "Phase B diagnostics require IDAES 2.12")
class PhaseBDiagnosticsTests(unittest.TestCase):
    def test_report_capture_preserves_successful_text(self):
        def callback(stream: io.StringIO):
            stream.write("diagnostic text")

        self.assertEqual(
            _capture_report(callback),
            {"status": "complete", "text": "diagnostic text"},
        )

    def test_report_capture_preserves_partial_failure(self):
        def callback(stream: io.StringIO):
            stream.write("partial")
            raise RuntimeError("missing backend")

        result = _capture_report(callback)
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["error_type"], "RuntimeError")
        self.assertEqual(result["text"], "partial")


if __name__ == "__main__":
    unittest.main()
