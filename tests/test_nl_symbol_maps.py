import hashlib
import json
import sys
import tempfile
import unittest
import weakref
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.compare_nl_symbol_maps import compare_sequences  # noqa: E402
from reproducibility.export_m5_nl_symbol_map import (  # noqa: E402
    _ordered_symbol_names,
    _inspect_nl,
    _resolve_symbol_object,
)
from reproducibility.run_m5_bakken_tax_series import _load_column_order  # noqa: E402


class _Named:
    def __init__(self, name):
        self.name = name


class _SymbolMap:
    def __init__(self, mapping):
        self.bySymbol = mapping


class _Model:
    def __init__(self, mapping):
        self.mapping = mapping

    def find_component(self, name):
        return self.mapping.get(name)


class NlSymbolMapTests(unittest.TestCase):
    def test_resolve_legacy_weak_reference(self):
        item = _Named("x")
        self.assertIs(_resolve_symbol_object(weakref.ref(item)), item)

    def test_ordered_symbol_names_uses_numeric_indices(self):
        symbol_map = _SymbolMap(
            {"v10": _Named("last"), "v2": _Named("middle"), "v0": _Named("first")}
        )
        result = _ordered_symbol_names(symbol_map)
        self.assertEqual(result["variables"], ["first", "middle", "last"])

    def test_compare_sequences_reports_displacement(self):
        result = compare_sequences(["a", "b", "c"], ["a", "c", "b"], top=2)
        self.assertTrue(result["component_sets_identical"])
        self.assertEqual(result["longest_common_prefix"], 1)
        self.assertEqual(result["same_position_count"], 1)
        self.assertEqual(result["maximum_absolute_position_displacement"], 1)

    def test_inspect_nl_streams_hash_and_header(self):
        content = b"first\nsecond\nthird\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.nl"
            path.write_bytes(content)
            result = _inspect_nl(path, header_lines=2)
        self.assertEqual(result["bytes"], len(content))
        self.assertEqual(result["sha256"], hashlib.sha256(content).hexdigest())
        self.assertEqual(result["header"], ["first", "second"])

    def test_load_column_order_resolves_and_verifies_names(self):
        names = ["model.x[1]", "model.x[0]"]
        digest = hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()
        first = _Named(names[0])
        second = _Named(names[1])
        payload = {
            "ordering": {"variables": {"names": names, "sha256": digest}}
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "map.json"
            path.write_text(json.dumps(payload))
            components, actual_digest = _load_column_order(
                _Model({names[0]: first, names[1]: second}), path
            )
        self.assertEqual(components, [first, second])
        self.assertEqual(actual_digest, digest)


if __name__ == "__main__":
    unittest.main()
