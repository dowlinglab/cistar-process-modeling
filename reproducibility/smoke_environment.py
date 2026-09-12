#!/usr/bin/env python3
"""Exercise the historical IDAES imports, M5 build, and HSL linear solvers."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import sys
from pathlib import Path

import idaes
import pandas as pd
import pyomo
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import Constraint, ConcreteModel, Objective, SolverFactory, Var, value


REPO_ROOT = Path(__file__).resolve().parents[1]
# Direct execution sets sys.path[0] to reproducibility/, not the repository
# root that contains the src package. Keep the documented command independent
# of the caller's PYTHONPATH.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

REPOSITORY_MODULES = (
    "src.dehydro_reactions",
    "src.rate_constant_custom_oligomerization",
    "src.rate_forms_oligo",
    "src.reaction_network_generator",
    "src.state_properties.properties_vap",
    "src.state_properties.properties_VLE_FpcTP",
    "src.state_properties.properties_vap_post_flash_2",
    "src.state_properties.properties_vap_H2_permeate",
    "src.unit_initialization",
    "src.costing_function",
    "src.emissions_calculations",
    "src.utility_minimization_1d",
)


def import_repository_modules() -> list[str]:
    imported = []
    for module_name in REPOSITORY_MODULES:
        importlib.import_module(module_name)
        imported.append(module_name)
    return imported


def build_m5_bakken() -> dict[str, int]:
    from src.unit_initialization import (
        create_flowsheet,
        define_arcs,
        define_models,
        set_scaling_factors,
        set_unit_model_variables,
    )

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
    return {
        "component_data_objects": len(list(model.component_data_objects())),
        "degrees_of_freedom": degrees_of_freedom(model),
    }


def check_hsl_solvers(ipopt: Path) -> dict[str, dict[str, float | str]]:
    # IDAES configures its own binary directory for subprocesses. On macOS that
    # can make an explicitly selected IPOPT executable load IDAES's different
    # libipopt ABI. Prefer the selected executable's sibling lib directory.
    solver_lib = str(ipopt.parent.parent / "lib")
    for variable in (
        "DYLD_LIBRARY_PATH",
        "DYLD_FALLBACK_LIBRARY_PATH",
        "LD_LIBRARY_PATH",
    ):
        current = os.environ.get(variable, "").strip(":")
        os.environ[variable] = f"{solver_lib}:{current}" if current else solver_lib
    results = {}
    for linear_solver in ("ma27", "ma57"):
        model = ConcreteModel()
        model.x = Var(initialize=0.5)
        model.y = Var(initialize=0.5)
        model.objective = Objective(
            expr=(model.x - 1) ** 2 + (model.y - 2) ** 2
        )
        model.constraint = Constraint(expr=model.x * model.y >= 1)
        solver = SolverFactory("ipopt", executable=str(ipopt))
        solver.options["linear_solver"] = linear_solver
        solve_result = solver.solve(model)
        results[linear_solver] = {
            "termination_condition": str(solve_result.solver.termination_condition),
            "objective": value(model.objective),
        }
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ipopt",
        type=Path,
        help="Optional HSL-enabled IPOPT executable to test with MA27 and MA57.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Only run imports/version checks (and the solver check, if requested).",
    )
    args = parser.parse_args()

    report = {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "idaes": idaes.__version__,
        "pyomo": pyomo.__version__,
        "imported_modules": import_repository_modules(),
    }
    if not args.skip_build:
        report["m5_bakken_build"] = build_m5_bakken()
    if args.ipopt:
        report["hsl_solver_smoke"] = check_hsl_solvers(args.ipopt.resolve())

    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
