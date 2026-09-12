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
    IDAES_CANDIDATE_COMMIT,
    _build_preoptimization_model,
    _collect_results,
    _configure_solver_environment,
    _solver_version,
    _write_report,
)
from run_m5_region import (
    PERTURBED_ZONES,
    PUBLISHED_TAX_USD_PER_KG,
    _base_optimum_region,
    _change_region,
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
        solve_result = solver.solve(model, tee=args.tee)
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
