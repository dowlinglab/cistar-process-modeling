#!/usr/bin/env python3
"""Audit a historical M5 regional pre-solve state without running IPOPT."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path
from typing import Any, Iterable

import idaes
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_diagnostics import (
    large_residuals_set,
    variables_near_bounds_set,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    badly_scaled_var_generator,
    unscaled_constraints_generator,
    unscaled_variables_generator,
)
from pyomo.environ import value

from run_m5_bakken_tax_series import (
    IDAES_CANDIDATE_COMMIT,
    _build_preoptimization_model,
    _checkpoint,
    _write_report,
)
from run_m5_region import (
    PUBLISHED_TAX_USD_PER_KG,
    REGIONS,
    _base_optimum_region,
    _change_region,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


def _component_record(component: Any) -> dict[str, Any]:
    return {"name": component.name, "value": value(component, exception=False)}


def _largest_residuals(model: Any, tolerance: float, limit: int) -> list[dict[str, Any]]:
    residuals = large_residuals_set(
        model, tol=tolerance, return_residual_values=True
    )
    ranked = sorted(residuals.items(), key=lambda item: item[1], reverse=True)
    return [
        {"name": constraint.name, "absolute_residual": residual}
        for constraint, residual in ranked[:limit]
    ]


def _badly_scaled_variables(
    model: Any, limit: int | None
) -> list[dict[str, Any]]:
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
    return records if limit is None else records[:limit]


def _near_bound_variables(model: Any, limit: int) -> list[dict[str, Any]]:
    variables = variables_near_bounds_set(model)
    records = [
        {
            "name": variable.name,
            "value": value(variable, exception=False),
            "lower_bound": value(variable.lb, exception=False),
            "upper_bound": value(variable.ub, exception=False),
        }
        for variable in variables
    ]
    return sorted(records, key=lambda record: record["name"])[:limit]


def _names(components: Iterable[Any], limit: int) -> tuple[int, list[str]]:
    names = sorted(component.name for component in components)
    return len(names), names[:limit]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", choices=REGIONS, required=True)
    parser.add_argument(
        "--initial-optimum",
        action="store_true",
        help="Audit the archived target-region optimum instead of the exact path.",
    )
    parser.add_argument("--residual-tolerance", type=float, default=1e-5)
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    model = _build_preoptimization_model(
        model_code=5,
        region="EF-Basin",
        costing_tax=PUBLISHED_TAX_USD_PER_KG,
        unit_initialization_region="Bakken",
    )
    initialization_chain = [
        "CISTAR_unit_initialization_Bakken_M5.json.gz",
        "CISTAR_solve_constrained_EF-Basin_M5_purge_0.01.json.gz",
        "CISTAR_solve_with_costing_EF-Basin_C_tax_0.045_M5_purge_0.01.json.gz",
    ]
    if args.region != "EF-Basin":
        initialization_chain.append(_change_region(model, args.region))

    if args.initial_optimum:
        checkpoint = _checkpoint(
            "CISTAR_optimal_solution_{}_C_tax_0.045_M5_purge_0.01_"
            "sequential_solve.json.gz".format(args.region)
        )
        ms.from_json(model, fname=str(checkpoint))
        initialization_chain.append(checkpoint.name)

    unfix_DOFs_pre_optimization(model)
    unscaled_variable_count, unscaled_variable_names = _names(
        unscaled_variables_generator(model), args.limit
    )
    unscaled_constraint_count, unscaled_constraint_names = _names(
        unscaled_constraints_generator(model), args.limit
    )
    near_bounds = variables_near_bounds_set(model)
    near_bound_count = len(near_bounds)
    badly_scaled = _badly_scaled_variables(model, limit=None)
    report = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "idaes_candidate_commit": IDAES_CANDIDATE_COMMIT,
            "pyomo": pyomo.__version__,
        },
        "case": {
            "model_code": 5,
            "region": args.region,
            "co2_tax_usd_per_kg": PUBLISHED_TAX_USD_PER_KG,
            "state": "archived_target_optimum"
            if args.initial_optimum
            else "exact_pre_solve_path",
            "base_optimum_region": _base_optimum_region(args.region),
        },
        "initialization_chain": initialization_chain,
        "degrees_of_freedom": degrees_of_freedom(model),
        "residual_tolerance": args.residual_tolerance,
        "largest_constraint_residuals": _largest_residuals(
            model, args.residual_tolerance, args.limit
        ),
        "badly_scaled_variable_count": len(badly_scaled),
        "badly_scaled_variables": badly_scaled[: args.limit],
        "near_bound_variable_count": near_bound_count,
        "near_bound_variables": _near_bound_variables(model, args.limit),
        "unscaled_variable_count": unscaled_variable_count,
        "unscaled_variables": unscaled_variable_names,
        "unscaled_constraint_count": unscaled_constraint_count,
        "unscaled_constraints": unscaled_constraint_names,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
