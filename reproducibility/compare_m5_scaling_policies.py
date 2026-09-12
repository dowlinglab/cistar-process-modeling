#!/usr/bin/env python3
"""Compare active NLP factors produced by modern and hybrid scaling policies."""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from scipy import sparse
from idaes.core.util.scaling import get_jacobian
from pyomo.environ import ConcreteModel, Constraint, Suffix, Var

from reproducibility.run_m5_bakken_tax_series import (
    _apply_active_nlp_autoscaling,
    _apply_modern_autoscaling,
    _build_preoptimization_model,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


def _capture_scaling_suffixes(model: Any) -> list[tuple[Any, list[tuple[Any, float]]]]:
    return [
        (suffix, list(suffix.items()))
        for suffix in model.component_objects(Suffix, descend_into=True)
        if suffix.local_name == "scaling_factor"
    ]


def _restore_scaling_suffixes(
    model: Any, state: list[tuple[Any, list[tuple[Any, float]]]]
) -> None:
    for suffix in model.component_objects(Suffix, descend_into=True):
        if suffix.local_name == "scaling_factor":
            suffix.clear()
    for suffix, entries in state:
        suffix.update(entries)


def _factor_maps(nlp: Any) -> dict[str, dict[str, float | None]]:
    from idaes.core.scaling.util import get_scaling_factor

    result: dict[str, dict[str, float | None]] = {
        "variables": {},
        "equalities": {},
        "inequalities": {},
    }
    for variable in nlp.vlist:
        factor = get_scaling_factor(variable)
        result["variables"][variable.name] = (
            None if factor is None else float(factor)
        )
    for constraint in nlp.clist:
        factor = get_scaling_factor(constraint)
        category = "equalities" if constraint.equality else "inequalities"
        result[category][constraint.name] = None if factor is None else float(factor)
    return result


def _minimal_autoscaler_probe() -> dict[str, Any]:
    """Expose equality/full-Jacobian row-index behavior on a three-row NLP."""
    from idaes.core.scaling import AutoScaler
    from idaes.core.scaling.util import get_scaling_factor

    model = ConcreteModel()
    model.x = Var(initialize=1.0)
    model.y = Var(initialize=1.0)
    model.eq1 = Constraint(expr=model.x == 1.0)
    model.ineq = Constraint(expr=10.0 * model.x <= 20.0)
    model.eq2 = Constraint(expr=100.0 * model.y == 100.0)
    jacobian, nlp = get_jacobian(model, scaled=False)
    AutoScaler(overwrite=True).scale_constraints_by_jacobian_norm(model, norm=2)
    return {
        "full_jacobian_rows": [
            {
                "index": index,
                "constraint": constraint.name,
                "equality": constraint.equality,
                "row_2_norm": float(sparse_norm),
            }
            for index, (constraint, sparse_norm) in enumerate(
                zip(
                    nlp.clist,
                    sparse.linalg.norm(jacobian, ord=2, axis=1),
                )
            )
        ],
        "eq2_expected_factor_from_own_row": 0.01,
        "eq2_observed_factor": float(get_scaling_factor(model.eq2)),
        "eq2_factor_corresponding_to_preceding_inequality_row": 0.1,
    }


def compare_factor_maps(
    reference: dict[str, float | None],
    candidate: dict[str, float | None],
    top: int = 20,
) -> dict[str, Any]:
    names = sorted(reference.keys() | candidate.keys())
    differences = []
    for name in names:
        first = reference.get(name)
        second = candidate.get(name)
        if first is None or second is None:
            log10_ratio = math.inf if first != second else 0.0
        else:
            log10_ratio = abs(math.log10(second / first))
        differences.append((log10_ratio, name, first, second))
    differences.sort(reverse=True)
    return {
        "reference_components": len(reference),
        "candidate_components": len(candidate),
        "component_names_identical": set(reference) == set(candidate),
        "reference_missing_factors": sum(value is None for value in reference.values()),
        "candidate_missing_factors": sum(value is None for value in candidate.values()),
        "counts_above_absolute_log10_ratio": {
            str(tolerance): sum(item[0] > tolerance for item in differences)
            for tolerance in (1e-12, 1e-9, 1e-6, 1e-3, 1.0)
        },
        "largest_factor_ratios": [
            {
                "name": name,
                "reference": first,
                "candidate": second,
                "absolute_log10_ratio": log10_ratio,
            }
            for log10_ratio, name, first, second in differences[:top]
        ],
    }


def compare_policies(norm: int = 2, top: int = 20) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    unfix_DOFs_pre_optimization(model)
    _, nlp = get_jacobian(model, scaled=False)
    inequality_positions = [
        {"index": index, "name": constraint.name}
        for index, constraint in enumerate(nlp.clist)
        if not constraint.equality
    ]
    inherited_state = _capture_scaling_suffixes(model)

    _apply_modern_autoscaling(model, norm)
    modern = _factor_maps(nlp)

    _restore_scaling_suffixes(model, inherited_state)
    hybrid_summary = _apply_active_nlp_autoscaling(
        model,
        norm,
        preserve_inherited_inequality_scaling=True,
    )
    hybrid = _factor_maps(nlp)

    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-SCALING-FACTOR-COMPARISON-001",
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "pyomo": pyomo.__version__,
        },
        "case": {"model_code": 5, "region": "Bakken", "co2_tax_usd_per_kg": 0.0},
        "reference_policy": "whole-model IDAES AutoScaler",
        "candidate_policy": "active-NLP rebuild preserving inherited inequalities",
        "hybrid_summary": hybrid_summary,
        "constraint_ordering": {
            "first_inequality_positions": inequality_positions[:10],
            "inequality_count": len(inequality_positions),
        },
        "minimal_autoscaler_probe": _minimal_autoscaler_probe(),
        "factor_comparisons": {
            category: compare_factor_maps(modern[category], hybrid[category], top)
            for category in ("variables", "equalities", "inequalities")
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--norm", type=int, default=2)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = compare_policies(arguments.norm, arguments.top)
    print(json.dumps(report, indent=2, sort_keys=True))
    _write_report(report, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
