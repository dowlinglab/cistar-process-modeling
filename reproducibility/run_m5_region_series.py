#!/usr/bin/env python3
"""Rerun the published M5 EF-zone series with incremental run records."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import idaes
import pyomo
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import SolverFactory

from run_m5_bakken_tax_series import (
    FileDeterminism,
    IDAES_CANDIDATE_COMMIT,
    _build_preoptimization_model,
    _collect_results,
    _configure_solver_environment,
    _load_column_order,
    _solver_version,
    _write_report,
)
from run_m5_region import (
    PERTURBED_ZONES,
    PUBLISHED_TAX_USD_PER_KG,
    _base_optimum_region,
    _change_region,
    _collect_figure_data,
    _load_archived_row,
)
from src.unit_initialization import (
    fix_DOFs_post_optimization,
    unfix_DOFs_pre_optimization,
)


ZONE_REGIONS = tuple(f"EF-{index}" for index in range(1, 13))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regions",
        nargs="+",
        choices=ZONE_REGIONS,
        default=list(ZONE_REGIONS),
    )
    parser.add_argument("--ipopt", type=Path, required=True)
    parser.add_argument(
        "--linear-solver", choices=("ma27", "ma57"), default="ma27"
    )
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--tee", action="store_true")
    parser.add_argument(
        "--inline-defined-variables",
        action="store_true",
        help="Set export_defined_variables=false for the modern Pyomo NL writer.",
    )
    parser.add_argument(
        "--file-determinism",
        choices=("ordered", "sort-indices", "sort-symbols"),
        default="ordered",
    )
    parser.add_argument("--column-order-from-symbol-map", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    started = time.time()
    ipopt = args.ipopt.resolve()
    solver_environment = _configure_solver_environment(ipopt)
    report = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "idaes_candidate_commit": IDAES_CANDIDATE_COMMIT,
            "pyomo": pyomo.__version__,
            "ipopt": _solver_version(ipopt, solver_environment),
            "linear_solver": args.linear_solver,
            "inline_defined_variables": args.inline_defined_variables,
            "file_determinism": args.file_determinism,
            "column_order_from_symbol_map": (
                str(args.column_order_from_symbol_map.resolve())
                if args.column_order_from_symbol_map is not None
                else None
            ),
        },
        "case": {
            "model_code": 5,
            "regions": args.regions,
            "co2_tax_usd_per_kg": PUBLISHED_TAX_USD_PER_KG,
        },
        "initialization_chain": [
            "CISTAR_unit_initialization_Bakken_M5.json.gz",
            "CISTAR_solve_constrained_EF-Basin_M5_purge_0.01.json.gz",
            "CISTAR_solve_with_costing_EF-Basin_C_tax_0.045_M5_purge_0.01.json.gz",
        ],
        "restart_rule": (
            "Each zone loads the archived EF-Basin optimum immediately before "
            "updating region data, except EF-9, which loads EF-8."
        ),
        "status": "building",
        "runs": [],
    }
    model = _build_preoptimization_model(
        model_code=5,
        region="EF-Basin",
        costing_tax=PUBLISHED_TAX_USD_PER_KG,
        unit_initialization_region="Bakken",
    )
    report["build_seconds"] = time.time() - started
    column_order = None
    if args.column_order_from_symbol_map is not None:
        column_order, column_digest = _load_column_order(
            model, args.column_order_from_symbol_map
        )
        report["environment"]["requested_column_order_sha256"] = column_digest
    report["status"] = "running"
    report["total_wall_seconds"] = time.time() - started
    _write_report(report, args.output)

    for region in args.regions:
        case_started = time.time()
        base_checkpoint = _change_region(model, region)
        unfix_DOFs_pre_optimization(model)
        perturbation = 653.0 if region in PERTURBED_ZONES else None
        if perturbation is not None:
            model.fs.H103.outlet.temperature.fix(perturbation)
            model.fs.H103.outlet.temperature.unfix()
        initial_dof = degrees_of_freedom(model)

        solver = SolverFactory("ipopt", executable=str(ipopt))
        solver.options.update(
            {
                "tol": 1e-6,
                "bound_push": 1e-8,
                "max_iter": args.max_iter,
                "linear_solver": args.linear_solver,
            }
        )
        try:
            writer_options = (
                {"export_defined_variables": False}
                if args.inline_defined_variables
                else {}
            )
            writer_options["file_determinism"] = {
                "ordered": FileDeterminism.ORDERED,
                "sort-indices": FileDeterminism.SORT_INDICES,
                "sort-symbols": FileDeterminism.SORT_SYMBOLS,
            }[args.file_determinism]
            if column_order is not None:
                writer_options["column_order"] = column_order
            solve_result = solver.solve(model, tee=args.tee, **writer_options)
        except KeyboardInterrupt:
            report["runs"].append(
                {
                    "region": region,
                    "base_optimum_region": _base_optimum_region(region),
                    "base_checkpoint": base_checkpoint,
                    "initial_temperature_perturbation_k": perturbation,
                    "initial_degrees_of_freedom": initial_dof,
                    "termination_condition": "interrupted_by_operator",
                    "wall_seconds": time.time() - case_started,
                }
            )
            report["status"] = "interrupted"
            report["total_wall_seconds"] = time.time() - started
            _write_report(report, args.output)
            raise
        except Exception as error:
            report["runs"].append(
                {
                    "region": region,
                    "base_optimum_region": _base_optimum_region(region),
                    "base_checkpoint": base_checkpoint,
                    "initial_temperature_perturbation_k": perturbation,
                    "initial_degrees_of_freedom": initial_dof,
                    "termination_condition": "python_exception",
                    "wall_seconds": time.time() - case_started,
                    "exception_type": type(error).__name__,
                    "exception_message": str(error),
                }
            )
            report["status"] = "error"
            report["total_wall_seconds"] = time.time() - started
            _write_report(report, args.output)
            raise
        fresh = _collect_results(model)
        archived = _load_archived_row(region)
        report["runs"].append(
            {
                "region": region,
                "base_optimum_region": _base_optimum_region(region),
                "base_checkpoint": base_checkpoint,
                "initial_temperature_perturbation_k": perturbation,
                "initial_degrees_of_freedom": initial_dof,
                "termination_condition": str(
                    solve_result.solver.termination_condition
                ),
                "wall_seconds": time.time() - case_started,
                "results_in_migrated_csv_units": fresh,
                "difference_from_migrated_csv": {
                    metric: fresh[metric] - archived[metric] for metric in fresh
                },
                "figure_data": _collect_figure_data(model),
            }
        )
        fix_DOFs_post_optimization(model)
        report["total_wall_seconds"] = time.time() - started
        _write_report(report, args.output)

    report["status"] = "complete"
    report["total_wall_seconds"] = time.time() - started
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0 if all(
        run["termination_condition"] == "optimal" for run in report["runs"]
    ) else 2


if __name__ == "__main__":
    raise SystemExit(main())
