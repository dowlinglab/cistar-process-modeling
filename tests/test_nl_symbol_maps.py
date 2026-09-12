import hashlib
import json
import sys
import tempfile
import unittest
import weakref
from pathlib import Path

from pyomo.environ import Block, ConcreteModel, Constraint, Set, Var


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from reproducibility.compare_nl_symbol_maps import compare_sequences  # noqa: E402
from reproducibility.compare_named_states import (  # noqa: E402
    _component_family,
    _unit_name,
    compare_values,
)
from reproducibility.export_m5_nl_symbol_map import (  # noqa: E402
    _ordered_symbol_names,
    _inspect_nl,
    _resolve_symbol_object,
)
from reproducibility.run_m5_bakken_tax_series import (  # noqa: E402
    _load_column_order,
    _regularize_pseudo_zero_inlet_phases,
)


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

    def test_compare_values_uses_bounded_scale(self):
        result = compare_values(100.0, 101.0)
        self.assertEqual(result["absolute_difference"], 1.0)
        self.assertAlmostEqual(result["scale_aware_difference"], 1.0 / 101.0)
        self.assertEqual(compare_values(0.0, 0.5)["scale_aware_difference"], 0.5)

    def test_state_aggregation_names(self):
        name = "fs.H106.properties[0.0].mole_frac_phase_comp[Vap,octane]"
        self.assertEqual(
            _component_family(name),
            "fs.H106.properties.mole_frac_phase_comp",
        )
        self.assertEqual(_unit_name(name), "fs.H106")

    def test_regularize_pseudo_zero_inlet_phases(self):
        model = ConcreteModel()
        model.fs = Block()
        for unit_name, phase, count, translator in (
            ("H105", "Liq", 2, False),
            ("H106", "Vap", 3, False),
            ("T103", "Liq", 2, True),
            ("T104", "Liq", 2, True),
        ):
            unit = Block()
            setattr(model.fs, unit_name, unit)
            if translator:
                unit.properties_in = Block([0.0])
                state = unit.properties_in[0.0]
            else:
                unit.control_volume = Block()
                unit.control_volume.properties_in = Block([0.0])
                state = unit.control_volume.properties_in[0.0]
            state.components = Set(initialize=[f"c{i}" for i in range(count)])
            state.phase_component_set = Set(
                dimen=2, initialize=[(phase, item) for item in state.components]
            )
            state.flow_mol_phase_comp = Var(
                state.phase_component_set, initialize=1e-8
            )
            state.mole_frac_phase_comp = Var(
                state.phase_component_set, initialize=1.0 / count
            )
            state.mole_frac_phase_comp_eq = Constraint(
                state.phase_component_set,
                rule=lambda block, p, j: block.flow_mol_phase_comp[p, j]
                == sum(
                    block.flow_mol_phase_comp[p, k] for k in block.components
                )
                * block.mole_frac_phase_comp[p, j],
            )
        result = _regularize_pseudo_zero_inlet_phases(model)
        self.assertEqual(result["states_regularized"], 4)
        self.assertEqual(result["constraints_deactivated"], 9)
        self.assertEqual(result["compositions_fixed"], 9)
        for detail in result["details"]:
            self.assertAlmostEqual(detail["fixed_fraction_sum"], 1.0)

    def test_tax_runner_exposes_opt_in_acceptable_termination(self):
        source = (
            REPO_ROOT / "reproducibility" / "run_m5_bakken_tax_series.py"
        ).read_text()
        self.assertIn('"--acceptable-tol"', source)
        self.assertIn('"--acceptable-iter"', source)
        self.assertIn('"--continue-after-failure"', source)


if __name__ == "__main__":
    unittest.main()
