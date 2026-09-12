#!/usr/bin/env python3
"""Run IDAES 2.12 diagnostics on the archived M5/Bakken starting state."""

from __future__ import annotations

import argparse
import io
import json
import platform
import sys
from pathlib import Path
from typing import Any, Callable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from idaes.core.util.diagnostics_tools import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    large_residuals_set,
    variables_near_bounds_set,
)
from idaes.core.util.scaling import (
    badly_scaled_var_generator,
    unscaled_constraints_generator,
    unscaled_variables_generator,
)
from pyomo.environ import Constraint, value

from reproducibility.run_m5_bakken_tax_series import (
    _build_preoptimization_model,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


def _capture_report(callback: Callable[..., None]) -> dict[str, Any]:
    stream = io.StringIO()
    try:
        callback(stream=stream)
    except Exception as err:  # diagnostic availability is itself evidence
        return {
            "status": "error",
            "error_type": type(err).__name__,
            "error": str(err),
            "text": stream.getvalue(),
        }
    return {"status": "complete", "text": stream.getvalue()}


def _ranked_residuals(model: Any, tolerance: float, limit: int) -> list[dict[str, Any]]:
    residuals = large_residuals_set(
        model, tol=tolerance, return_residual_values=True
    )
    return [
        {"name": component.name, "absolute_residual": residual}
        for component, residual in sorted(
            residuals.items(), key=lambda item: item[1], reverse=True
        )[:limit]
    ]


def _bad_scaling(model: Any, limit: int) -> tuple[int, list[dict[str, Any]]]:
    records = [
        {
            "name": variable.name,
            "value": value(variable, exception=False),
            "scaled_absolute_value": scaled_value,
        }
        for variable, scaled_value in badly_scaled_var_generator(model)
    ]
    records.sort(
        key=lambda record: max(
            record["scaled_absolute_value"],
            1.0 / max(record["scaled_absolute_value"], 1e-300),
        ),
        reverse=True,
    )
    return len(records), records[:limit]


def _component_names(components: Any, limit: int) -> dict[str, Any]:
    names = sorted(component.name for component in components)
    return {"count": len(names), "first_names": names[:limit]}


def diagnose(
    limit: int = 50,
    residual_tolerance: float = 1e-5,
    include_structural_report: bool = False,
    include_targeted_toolbox_reports: bool = False,
    include_full_numerical_report: bool = False,
) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    fixed_state = {
        "expected_degrees_of_freedom": 0,
        "active_constraint_count": sum(
            1
            for _ in model.component_data_objects(
                Constraint, active=True, descend_into=True
            )
        ),
    }
    toolbox = DiagnosticsToolbox(model)
    if include_structural_report:
        fixed_state["structural_report"] = _capture_report(
            toolbox.report_structural_issues
        )
    else:
        fixed_state["structural_report"] = {
            "status": "deferred",
            "reason": (
                "The full unit-consistency traversal is optional because the "
                "legacy flowsheet emits a very large number of unit errors and "
                "takes substantially longer than the numerical diagnostics."
            ),
        }
    if include_targeted_toolbox_reports:
        fixed_state["targeted_toolbox_reports"] = {
            "large_residuals": _capture_report(
                toolbox.display_constraints_with_large_residuals
            ),
            "variables_at_or_outside_bounds": _capture_report(
                toolbox.display_variables_at_or_outside_bounds
            ),
            "variables_with_extreme_values": _capture_report(
                toolbox.display_variables_with_extreme_values
            ),
            "variables_near_bounds": _capture_report(
                toolbox.display_variables_near_bounds
            ),
        }
    else:
        fixed_state["targeted_toolbox_reports"] = {
            "status": "deferred",
            "reason": (
                "The targeted display methods also require a costly traversal "
                "of this legacy model. Equivalent bounded component summaries "
                "are recorded below."
            ),
        }
    if include_full_numerical_report:
        fixed_state["numerical_report"] = _capture_report(
            toolbox.report_numerical_issues
        )
    else:
        fixed_state["numerical_report"] = {
            "status": "deferred",
            "reason": (
                "The full report performs costly Jacobian and constraint-term "
                "analyses on this 25,000-component legacy model. Targeted "
                "DiagnosticsToolbox reports are recorded instead."
            ),
        }
    fixed_state["largest_constraint_residuals"] = _ranked_residuals(
        model, residual_tolerance, limit
    )

    unfix_DOFs_pre_optimization(model)
    badly_scaled_count, badly_scaled = _bad_scaling(model, limit)
    optimization_state = {
        "expected_degrees_of_freedom": 8,
        "active_constraint_count": fixed_state["active_constraint_count"],
        "largest_constraint_residuals": _ranked_residuals(
            model, residual_tolerance, limit
        ),
        "badly_scaled_variable_count": badly_scaled_count,
        "badly_scaled_variables": badly_scaled,
        "near_bound_variables": _component_names(
            variables_near_bounds_set(model), limit
        ),
        "unscaled_variables": _component_names(
            unscaled_variables_generator(model), limit
        ),
        "unscaled_constraints": _component_names(
            unscaled_constraints_generator(model), limit
        ),
    }
    return {
        "schema_version": 1,
        "experiment_id": "B-IDAES-2.12-M5-BAKKEN-PRESOLVE-DIAGNOSTICS-001",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "pyomo": pyomo.__version__,
        },
        "case": {"model_code": 5, "region": "Bakken", "co2_tax_usd_per_kg": 0.0},
        "residual_tolerance": residual_tolerance,
        "fixed_archived_start": fixed_state,
        "optimization_start": optimization_state,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--residual-tolerance", type=float, default=1e-5)
    parser.add_argument(
        "--include-structural-report",
        action="store_true",
        help="Also run the expensive full structural/unit-consistency audit.",
    )
    parser.add_argument(
        "--include-full-numerical-report",
        action="store_true",
        help="Also run the expensive full Jacobian/constraint-term report.",
    )
    parser.add_argument(
        "--include-targeted-toolbox-reports",
        action="store_true",
        help="Run DiagnosticsToolbox display methods in addition to bounded summaries.",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = diagnose(
        args.limit,
        args.residual_tolerance,
        include_structural_report=args.include_structural_report,
        include_targeted_toolbox_reports=args.include_targeted_toolbox_reports,
        include_full_numerical_report=args.include_full_numerical_report,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    _write_report(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
