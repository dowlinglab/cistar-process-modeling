#!/usr/bin/env python3
"""Rerun one published M5 regional optimization without overwriting archives."""

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
import pandas as pd
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import arcs_to_stream_dict, create_stream_table_dataframe
from pyomo.environ import SolverFactory, value

from run_m5_bakken_tax_series import (
    IDAES_CANDIDATE_COMMIT,
    REPO_ROOT,
    _build_preoptimization_model,
    _checkpoint,
    _collect_results,
    _configure_solver_environment,
    _solver_version,
    _write_report,
)
from src.emissions_calculations import (
    calc_lhv_values,
    delete_region_specific_components,
)
from src.unit_initialization import (
    fix_DOFs_post_optimization,
    unfix_DOFs_pre_optimization,
)
from src.utility_minimization_1d import gen_curves, heat_ex_data, return_HX_results


PUBLISHED_TAX_USD_PER_KG = 0.045
REGIONS = ("EF-Basin",) + tuple(f"EF-{index}" for index in range(1, 13))
PERTURBED_ZONES = {"EF-1", "EF-6", "EF-7", "EF-10", "EF-11"}
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


def _load_archived_row(region: str) -> dict[str, float]:
    path = REPO_ROOT / "results" / "optimal_data_wrt_region.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["ROK_model"] == "M5" and row["Region"] == region:
                return {metric: float(row[metric]) for metric in RESULT_METRICS}
    raise KeyError(f"No archived M5/{region} result in {path}")


def _base_optimum_region(region: str) -> str:
    return "EF-8" if region == "EF-9" else "EF-Basin"


def _change_region(model: Any, region: str) -> str:
    base_region = _base_optimum_region(region)
    checkpoint = _checkpoint(
        "CISTAR_optimal_solution_{}_C_tax_0.045_M5_purge_0.01_"
        "sequential_solve.json.gz".format(base_region)
    )
    ms.from_json(model, fname=str(checkpoint))

    inlet_data = pd.read_csv(REPO_ROOT / "data" / "NGL_compositions.csv")
    for _, row in inlet_data.iterrows():
        variable = model.fs.M101.feed.mole_frac_comp[0, row["Species"]]
        variable.unfix()
        value = 1e-6 if row[region] == 0.0 else round(row[region], 4)
        variable.fix(value)

    delete_region_specific_components(model)
    calc_lhv_values(
        model,
        region,
        str(REPO_ROOT / "data" / "LHV.xlsx"),
        str(REPO_ROOT / "data" / "NGL_compositions.csv"),
        str(REPO_ROOT / "data" / "NGL_fraction.csv"),
    )
    emissions = pd.read_csv(REPO_ROOT / "data" / "emissions_factor_by_region.csv")
    model.fs.upstream_emission_factor = float(emissions[region].iloc[0])
    return checkpoint.name


def _base_report(
    region: str,
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
            "model_code": 5,
            "region": region,
            "co2_tax_usd_per_kg": PUBLISHED_TAX_USD_PER_KG,
        },
        "initialization_chain": [
            "CISTAR_unit_initialization_Bakken_M5.json.gz",
            "CISTAR_solve_constrained_EF-Basin_M5_purge_0.01.json.gz",
            "CISTAR_solve_with_costing_EF-Basin_C_tax_0.045_M5_purge_0.01.json.gz",
        ],
        "status": "building",
    }


def _dataframe_payload(frame: pd.DataFrame) -> dict[str, Any]:
    """Return a JSON-safe, orientation-preserving DataFrame representation."""
    return json.loads(frame.to_json(orient="split", double_precision=15))


def _collect_figure_data(model: Any) -> dict[str, Any]:
    """Collect the numerical series used by the regional analysis figures."""
    heating = [model.fs.H101, model.fs.H103, model.fs.R101]
    cooling = [
        model.fs.H102,
        model.fs.H104,
        model.fs.H105,
        model.fs.H106,
        model.fs.R102,
    ]
    curve_data = heat_ex_data(model.fs, heating, cooling)
    hot_temperature, hot_heat = gen_curves(
        curve_data.Cooling_Tin,
        curve_data.Cooling_Tout,
        curve_data.Cooling_Q,
    )
    cold_temperature, cold_heat = gen_curves(
        curve_data.Heating_Tin,
        curve_data.Heating_Tout,
        -curve_data.Heating_Q,
    )
    cold_heat = cold_heat + sum(curve_data.Heating_Q) + value(curve_data.Qw)

    component_flow = {
        component: value(
            model.fs.F102.liq_outlet.flow_mol_phase_comp[0, "Liq", component]
        )
        for component in model.fs.liquid
    }
    total_liquid_flow = sum(component_flow.values())

    return {
        "upstream_emissions_kg_co2e_per_gj": value(
            model.fs.upstream_emissions
        ),
        "liquid_product": {
            "total_flow_mol_per_s": total_liquid_flow,
            "component_flow_mol_per_s": component_flow,
            "component_mole_percent": {
                component: 100.0 * flow / total_liquid_flow
                for component, flow in component_flow.items()
            },
            "component_lhv_contribution_mj_per_s": {
                component: value(
                    model.fs.LHV_per_component_per_stream[
                        "liq_outlet", component
                    ]
                )
                * total_liquid_flow
                for component in model.fs.liquid
            },
        },
        "composite_curves": {
            "hot": {
                "temperature_k": hot_temperature.tolist(),
                "cumulative_heat_gj_per_hour": (-hot_heat).tolist(),
            },
            "cold": {
                "temperature_k": cold_temperature.tolist(),
                "cumulative_heat_gj_per_hour": cold_heat.tolist(),
            },
        },
        "heat_exchanger_table": _dataframe_payload(
            return_HX_results(model.fs, heating + cooling)
        ),
        "stream_table": _dataframe_payload(
            create_stream_table_dataframe(
                arcs_to_stream_dict(model, descend_into=True)
            )
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", choices=REGIONS, required=True)
    parser.add_argument("--ipopt", type=Path, required=True)
    parser.add_argument(
        "--linear-solver", choices=("ma27", "ma57"), default="ma27"
    )
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument(
        "--initial-optimum",
        action="store_true",
        help="Load the archived target-region optimum before solving.",
    )
    parser.add_argument("--tee", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    started = time.time()
    ipopt = args.ipopt.resolve()
    solver_environment = _configure_solver_environment(ipopt)
    report = _base_report(
        args.region, ipopt, args.linear_solver, solver_environment
    )
    model = _build_preoptimization_model(
        model_code=5,
        region="EF-Basin",
        costing_tax=PUBLISHED_TAX_USD_PER_KG,
        unit_initialization_region="Bakken",
    )

    if args.region == "EF-Basin":
        perturbation = 593.0
    else:
        base_checkpoint = _change_region(model, args.region)
        report["initialization_chain"].append(base_checkpoint)
        report["base_optimum_region"] = _base_optimum_region(args.region)
        perturbation = 653.0 if args.region in PERTURBED_ZONES else None

    if args.initial_optimum:
        checkpoint = _checkpoint(
            "CISTAR_optimal_solution_{}_C_tax_0.045_M5_purge_0.01_"
            "sequential_solve.json.gz".format(args.region)
        )
        ms.from_json(model, fname=str(checkpoint))
        report["initialization_chain"].append(checkpoint.name)
        report["archived_initial_optimum"] = checkpoint.name
        perturbation = None

    report["build_seconds"] = time.time() - started
    report["status"] = "running"
    report["total_wall_seconds"] = time.time() - started
    _write_report(report, args.output)

    unfix_DOFs_pre_optimization(model)
    if perturbation is not None:
        model.fs.H103.outlet.temperature.fix(perturbation)
        model.fs.H103.outlet.temperature.unfix()
        report["initial_temperature_perturbation"] = {
            "variable": "fs.H103.outlet.temperature",
            "value_k": perturbation,
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
    try:
        solve_result = solver.solve(model, tee=args.tee)
    except KeyboardInterrupt:
        report["run"] = {
            "initial_degrees_of_freedom": initial_dof,
            "termination_condition": "interrupted_by_operator",
            "wall_seconds": time.time() - solve_started,
        }
        report["status"] = "interrupted"
        report["total_wall_seconds"] = time.time() - started
        _write_report(report, args.output)
        raise
    except Exception as error:
        report["run"] = {
            "initial_degrees_of_freedom": initial_dof,
            "termination_condition": "python_exception",
            "wall_seconds": time.time() - solve_started,
            "exception_type": type(error).__name__,
            "exception_message": str(error),
        }
        report["status"] = "error"
        report["total_wall_seconds"] = time.time() - started
        _write_report(report, args.output)
        raise
    fresh = _collect_results(model)
    archived = _load_archived_row(args.region)
    report["run"] = {
        "initial_degrees_of_freedom": initial_dof,
        "termination_condition": str(solve_result.solver.termination_condition),
        "wall_seconds": time.time() - solve_started,
        "results_in_migrated_csv_units": fresh,
        "difference_from_migrated_csv": {
            metric: fresh[metric] - archived[metric] for metric in fresh
        },
        "figure_data": _collect_figure_data(model),
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
