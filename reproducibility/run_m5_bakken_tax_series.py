#!/usr/bin/env python3
"""Rerun the published M5/Bakken carbon-tax sequence without overwriting archives."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import idaes
import pandas as pd
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import SolverFactory, value

from src.costing_function import add_costing, calculate_costs_for_objective
from src.emissions_calculations import (
    calc_lhv_values,
    calculate_emissions,
    calculate_stream_energies,
    create_ghg_objective,
)
from src.unit_initialization import (
    create_flowsheet,
    define_arcs,
    define_models,
    fix_DOFs_post_optimization,
    set_scaling_factors,
    set_unit_model_variables,
    unfix_DOFs_pre_optimization,
    update_model_after_initialization,
    update_model_for_optimization,
    vapor_only_to_vapor_liquid_reformulate,
)
from src.utility_minimization_1d import min_utility


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TAX_RATES = (0.0, 1e-5, 1e-3, 1.7e-2, 4.5e-2, 1.9e-1, 4.1e-1)
IDAES_CANDIDATE_COMMIT = "66935c80a5aafc3ffc9ab3d387e488cddd4f233b"


def _checkpoint(name: str) -> Path:
    return REPO_ROOT / "initialization_files" / name


def _build_preoptimization_model() -> Any:
    model = create_flowsheet(5)
    define_models(model, catalyst_mass=1167.003367)
    define_arcs(model)

    inlet_data = pd.read_csv(REPO_ROOT / "data" / "NGL_compositions.csv")
    composition = {
        row["Species"]: (1e-6 if row["Bakken"] == 0.0 else round(row["Bakken"], 4))
        for _, row in inlet_data.iterrows()
    }
    conversion = {"ethane": 0.3566, "propane": 0.6632, "nbutane": 0.5188}
    set_unit_model_variables(
        model,
        model_code=5,
        feed_flow_rate=481.3888889,
        feed_temp=308.0,
        feed_pressure=700000.0,
        inlet_composition_dict=composition,
        dehydro_conv_dict=conversion,
    )
    set_scaling_factors(
        model,
        flow_mol_scaling_factor=1e-3,
        inlet_composition_dict=composition,
    )

    ms.from_json(
        model,
        fname=str(_checkpoint("CISTAR_unit_initialization_Bakken_M5.json.gz")),
    )
    update_model_after_initialization(model)
    # Preserve the published notebook sequence exactly. The second call is
    # idempotent and reports that no degenerate constraints remain.
    vapor_only_to_vapor_liquid_reformulate(model.fs.T102)
    vapor_only_to_vapor_liquid_reformulate(model.fs.T102)
    ms.from_json(
        model,
        fname=str(
            _checkpoint("CISTAR_solve_constrained_Bakken_M5_purge_0.01.json.gz")
        ),
    )

    add_costing(model)
    min_utility(
        model.fs,
        [model.fs.H101, model.fs.H103, model.fs.R101],
        [model.fs.H102, model.fs.H104, model.fs.H105, model.fs.H106, model.fs.R102],
        10.0,
    )
    model.fs.Qs.fix()
    calc_lhv_values(
        model,
        "Bakken",
        str(REPO_ROOT / "data" / "LHV.xlsx"),
        str(REPO_ROOT / "data" / "NGL_compositions.csv"),
        str(REPO_ROOT / "data" / "NGL_fraction.csv"),
    )
    calculate_stream_energies(model)
    calculate_emissions(
        model,
        "Bakken",
        str(REPO_ROOT / "data" / "emissions_factor_by_region.csv"),
    )
    create_ghg_objective(model)
    calculate_costs_for_objective(model, c_tax_flag=True, c_tax_val=0.0)
    ms.from_json(
        model,
        fname=str(
            _checkpoint(
                "CISTAR_solve_with_costing_Bakken_C_tax_0.0_M5_purge_0.01.json.gz"
            )
        ),
    )
    update_model_for_optimization(model)
    return model


def _configure_solver_environment(ipopt: Path) -> dict[str, str]:
    environment = os.environ.copy()
    solver_lib = str(ipopt.parent.parent / "lib")
    for variable in (
        "DYLD_LIBRARY_PATH",
        "DYLD_FALLBACK_LIBRARY_PATH",
        "LD_LIBRARY_PATH",
    ):
        current = environment.get(variable, "").strip(":")
        environment[variable] = (
            f"{solver_lib}:{current}" if current else solver_lib
        )
        os.environ[variable] = environment[variable]
    return environment


def _solver_version(ipopt: Path, environment: dict[str, str]) -> str:
    completed = subprocess.run(
        [str(ipopt), "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    return completed.stdout.strip()


def _collect_results(model: Any) -> dict[str, float]:
    return {
        "MSP": value(model.fs.min_sell_price) / 1000,
        "Downstream-em": value(model.fs.downstream_emissions),
        "Product-LHV": value(
            model.fs.fuel_energy + model.fs.gas_energy + model.fs.h2_energy
        )
        * 1000
        / 3600,
        "H2-rebate": value(model.fs.H2_purge_sell_price),
        "TAC": value(model.fs.TAC),
        "T_R102": value(model.fs.H103.outlet.temperature[0]),
        "T_H104": value(model.fs.H104.outlet.temperature[0]),
        "T_H105": value(model.fs.H105.outlet.temperature[0]),
        "P_F101": value(model.fs.F101.liq_outlet.pressure[0]),
        "T_H106": value(model.fs.H106.outlet.temperature[0]),
        "P_H106": value(model.fs.H106.outlet.pressure[0]),
        "P_F102": value(model.fs.F102.liq_outlet.pressure[0]),
        "Qs": value(model.fs.Qs) * 1000 / 3600,
        "Qw": value(model.fs.Qw) * 1000 / 3600,
    }


def _load_archived_rows() -> dict[float, dict[str, float]]:
    path = REPO_ROOT / "results" / "optimal_data_wrt_c_tax_rates.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        rows = csv.DictReader(stream)
        return {
            float(row["C-tax-rate"]): {
                metric: float(row[metric])
                for metric in (
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
            }
            for row in rows
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ipopt", type=Path, required=True)
    parser.add_argument(
        "--linear-solver", choices=("ma27", "ma57"), default="ma27"
    )
    parser.add_argument(
        "--tax-rates",
        type=float,
        nargs="+",
        default=list(DEFAULT_TAX_RATES),
        metavar="USD_PER_KG",
    )
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument(
        "--initial-optimal-tax",
        type=float,
        metavar="USD_PER_KG",
        help=(
            "Load the archived M5/Bakken optimum at this tax before the first "
            "requested solve. Useful for isolating one sequential transition."
        ),
    )
    parser.add_argument("--tee", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    started = time.time()
    ipopt = args.ipopt.resolve()
    solver_environment = _configure_solver_environment(ipopt)
    report: dict[str, Any] = {
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
        "case": {"model_code": 5, "region": "Bakken"},
        "tax_rate_units": "USD/kg CO2e",
        "initialization_chain": [
            "CISTAR_unit_initialization_Bakken_M5.json.gz",
            "CISTAR_solve_constrained_Bakken_M5_purge_0.01.json.gz",
            "CISTAR_solve_with_costing_Bakken_C_tax_0.0_M5_purge_0.01.json.gz",
        ],
        "runs": [],
    }

    model = _build_preoptimization_model()
    if args.initial_optimal_tax is not None:
        checkpoint = _checkpoint(
            "CISTAR_optimal_solution_Bakken_C_tax_{}_M5_purge_0.01_"
            "sequential_solve.json.gz".format(args.initial_optimal_tax)
        )
        if not checkpoint.exists():
            raise FileNotFoundError(f"Archived starting point not found: {checkpoint}")
        model.fs.c_tax_rate = args.initial_optimal_tax
        unfix_DOFs_pre_optimization(model)
        ms.from_json(model, fname=str(checkpoint))
        report["archived_initial_optimum"] = checkpoint.name
    report["build_seconds"] = time.time() - started
    archived_rows = _load_archived_rows()

    for tax_rate in args.tax_rates:
        case_started = time.time()
        model.fs.c_tax_rate = tax_rate
        unfix_DOFs_pre_optimization(model)
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
        archived = archived_rows.get(tax_rate * 1000)
        comparison = None
        if archived is not None:
            comparison = {
                metric: fresh[metric] - archived[metric] for metric in fresh
            }
        report["runs"].append(
            {
                "co2_tax_usd_per_kg": tax_rate,
                "initial_degrees_of_freedom": initial_dof,
                "termination_condition": str(
                    solve_result.solver.termination_condition
                ),
                "wall_seconds": time.time() - case_started,
                "results_in_migrated_csv_units": fresh,
                "difference_from_migrated_csv": comparison,
            }
        )
        fix_DOFs_post_optimization(model)

    report["total_wall_seconds"] = time.time() - started
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if all(
        run["termination_condition"] == "optimal" for run in report["runs"]
    ) else 2


if __name__ == "__main__":
    raise SystemExit(main())
