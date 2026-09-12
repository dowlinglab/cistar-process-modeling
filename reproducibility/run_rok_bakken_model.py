#!/usr/bin/env python3
"""Rerun one published M2-M4 Bakken optimization without overwriting archives."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import idaes
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import SolverFactory

from run_m5_bakken_tax_series import (
    IDAES_CANDIDATE_COMMIT,
    REPO_ROOT,
    _build_preoptimization_model,
    _collect_results,
    _configure_solver_environment,
    _solver_version,
    _write_report,
)
from src.unit_initialization import (
    fix_DOFs_post_optimization,
    unfix_DOFs_pre_optimization,
)
from src.result_extraction import collect_figure_data


PUBLISHED_TAX_USD_PER_KG = 0.045
RESULT_METRICS = (
    "MSP",
    "Downstream-em",
    "Product-LHV",
    "H2-rebate",
    "TAC",
    "T_R102",
    "T_H104",
    "T_H105",
    "P_F101",
    "T_H106",
    "P_H106",
    "P_F102",
    "Qs",
    "Qw",
)


def _load_archived_row(model_code: int) -> dict[str, float]:
    path = REPO_ROOT / "results" / "optimal_data_wrt_ROK_models.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["ROK_model"] == f"M{model_code}" and row["Region"] == "Bakken":
                return {metric: float(row[metric]) for metric in RESULT_METRICS}
    raise KeyError(f"No archived M{model_code}/Bakken result in {path}")


def _base_report(
    model_code: int,
    ipopt: Path,
    linear_solver: str,
    solver_environment: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "idaes_candidate_commit": IDAES_CANDIDATE_COMMIT,
            "pyomo": pyomo.__version__,
            "ipopt": _solver_version(ipopt, solver_environment),
            "linear_solver": linear_solver,
        },
        "case": {
            "model_code": model_code,
            "region": "Bakken",
            "co2_tax_usd_per_kg": PUBLISHED_TAX_USD_PER_KG,
        },
        "initialization_chain": [
            f"CISTAR_unit_initialization_Bakken_M{model_code}.json.gz",
            f"CISTAR_solve_constrained_Bakken_M{model_code}_purge_0.01.json.gz",
            (
                "CISTAR_solve_with_costing_Bakken_C_tax_0.045_"
                f"M{model_code}_purge_0.01.json.gz"
            ),
        ],
        "status": "building",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-code", type=int, choices=(2, 3, 4), required=True)
    parser.add_argument("--ipopt", type=Path, required=True)
    parser.add_argument(
        "--linear-solver", choices=("ma27", "ma57"), default="ma27"
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=500,
        help="IPOPT iteration cap; the published M2-M4 notebooks use 500.",
    )
    parser.add_argument(
        "--initial-optimum",
        action="store_true",
        help="Load the archived optimum after rebuilding the costed model.",
    )
    parser.add_argument("--tee", action="store_true")
    parser.add_argument(
        "--inline-defined-variables",
        action="store_true",
        help="Set export_defined_variables=false for the modern Pyomo NL writer.",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    started = time.time()
    ipopt = args.ipopt.resolve()
    solver_environment = _configure_solver_environment(ipopt)
    report = _base_report(
        args.model_code, ipopt, args.linear_solver, solver_environment
    )
    report["environment"]["inline_defined_variables"] = (
        args.inline_defined_variables
    )

    model = _build_preoptimization_model(
        model_code=args.model_code,
        region="Bakken",
        costing_tax=PUBLISHED_TAX_USD_PER_KG,
    )
    if args.initial_optimum:
        checkpoint = (
            REPO_ROOT
            / "initialization_files"
            / (
                "CISTAR_optimal_solution_Bakken_C_tax_0.045_"
                f"M{args.model_code}_purge_0.01.json.gz"
            )
        )
        ms.from_json(model, fname=str(checkpoint))
        report["archived_initial_optimum"] = checkpoint.name
    report["build_seconds"] = time.time() - started
    report["status"] = "running"
    report["total_wall_seconds"] = time.time() - started
    _write_report(report, args.output)

    unfix_DOFs_pre_optimization(model)
    # Preserve the notebook's explicit convergence perturbation. Fixing sets
    # the current value to 590 K; unfixing immediately afterward leaves the
    # variable free but retains that initial value for the NLP solve.
    model.fs.H103.outlet.temperature.fix(590.0)
    model.fs.H103.outlet.temperature.unfix()
    report["initial_temperature_perturbation"] = {
        "variable": "fs.H103.outlet.temperature",
        "value_k": 590.0,
        "method": "fix_then_unfix",
    }
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
    solve_started = time.time()
    writer_options = (
        {"export_defined_variables": False}
        if args.inline_defined_variables
        else {}
    )
    solve_result = solver.solve(model, tee=args.tee, **writer_options)
    fresh = _collect_results(model)
    archived = _load_archived_row(args.model_code)
    report["run"] = {
        "initial_degrees_of_freedom": initial_dof,
        "termination_condition": str(solve_result.solver.termination_condition),
        "wall_seconds": time.time() - solve_started,
        "results_in_migrated_csv_units": fresh,
        "difference_from_migrated_csv": {
            metric: fresh[metric] - archived[metric] for metric in fresh
        },
        "figure_data": collect_figure_data(model),
    }
    fix_DOFs_post_optimization(model)
    report["status"] = "complete"
    report["total_wall_seconds"] = time.time() - started
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0 if report["run"]["termination_condition"] == "optimal" else 2


if __name__ == "__main__":
    raise SystemExit(main())
