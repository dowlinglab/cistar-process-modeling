#!/usr/bin/env python3
"""Rerun the published M5/Bakken carbon-tax sequence without overwriting archives."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import numpy as np
import pandas as pd
import pyomo
from scipy import sparse
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    get_constraint_transform_applied_scaling_factor,
    get_jacobian,
    set_scaling_factor,
)
from pyomo.environ import Constraint, SolverFactory, Suffix, Var, value
try:
    from pyomo.repn.util import FileDeterminism
except ImportError:  # Pyomo 6.4 keeps this enum in the NL writer module
    from pyomo.repn.plugins.nl_writer import FileDeterminism

from src.costing_function import add_costing, calculate_costs_for_objective
from src.emissions_calculations import (
    calc_lhv_values,
    calculate_emissions,
    calculate_stream_energies,
    create_ghg_objective,
)
from src.result_extraction import collect_figure_data
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


DEFAULT_TAX_RATES = (0.0, 1e-5, 1e-3, 1.7e-2, 4.5e-2, 1.9e-1, 4.1e-1)
IDAES_CANDIDATE_COMMIT = "66935c80a5aafc3ffc9ab3d387e488cddd4f233b"
DESIGN_VARIABLE_NAMES = (
    "Qs",
    "H103_temperature",
    "H104_temperature",
    "H105_temperature",
    "F101_deltaP",
    "H106_temperature",
    "H106_pressure",
    "F102_deltaP",
)


def _idaes_provenance() -> dict[str, str]:
    """Describe the installed IDAES without mislabeling modern releases."""
    provenance = {"idaes": idaes.__version__}
    if idaes.__version__.startswith("2.0.0.dev3"):
        provenance["idaes_candidate_commit"] = IDAES_CANDIDATE_COMMIT
    else:
        provenance["historical_reference_commit"] = IDAES_CANDIDATE_COMMIT
    return provenance


def _apply_modern_autoscaling(model: Any, norm: int) -> None:
    """Overwrite legacy suffixes using the IDAES 2.12 scaling toolbox."""
    try:
        from idaes.core.scaling import AutoScaler
    except ImportError as err:
        raise RuntimeError(
            "Modern autoscaling requires an IDAES release that exposes AutoScaler."
        ) from err
    AutoScaler(overwrite=True).scale_model(model, norm=norm)


def _clear_scaling_suffixes(model: Any) -> dict[str, int]:
    """Clear inherited scaling metadata without changing transformed equations."""
    suffixes = 0
    entries = 0
    for suffix in model.component_objects(Suffix, descend_into=True):
        if suffix.local_name == "scaling_factor":
            suffixes += 1
            entries += len(suffix)
            suffix.clear()
    return {"suffixes_cleared": suffixes, "entries_cleared": entries}


def _apply_active_nlp_autoscaling(
    model: Any,
    norm: int,
    preserve_inherited_inequality_scaling: bool = False,
) -> dict[str, int]:
    """Rebuild suffix scaling only for variables/constraints in the active NLP."""
    try:
        from idaes.core.scaling import AutoScaler
        from idaes.core.scaling.util import set_scaling_factor
    except ImportError as err:
        raise RuntimeError(
            "Active-NLP autoscaling requires an IDAES release that exposes AutoScaler."
        ) from err

    inherited_inequality_factors: dict[str, float] = {}
    if preserve_inherited_inequality_scaling:
        from idaes.core.scaling.util import get_scaling_factor

        for constraint in model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            if not constraint.equality:
                factor = get_scaling_factor(constraint)
                inherited_inequality_factors[constraint.name] = (
                    1.0 if factor is None else float(factor)
                )

    cleared = _clear_scaling_suffixes(model)
    _, nlp = get_jacobian(model, scaled=False)
    active_variables = list(nlp.vlist)
    active_constraints = list(nlp.clist)

    scaler = AutoScaler(overwrite=True)
    for variable in active_variables:
        scaler.scale_variables_by_magnitude(variable)
    scaled_jacobian, scaled_nlp = get_jacobian(model, scaled=True)
    if [item.name for item in scaled_nlp.clist] != [
        item.name for item in active_constraints
    ]:
        raise RuntimeError("Active constraint ordering changed during scaling.")
    row_norms = np.asarray(
        sparse.linalg.norm(scaled_jacobian, ord=int(norm), axis=1)
    ).reshape(-1)
    for constraint, row_norm in zip(active_constraints, row_norms):
        if constraint.name in inherited_inequality_factors:
            factor = inherited_inequality_factors[constraint.name]
        else:
            factor = (
                1.0
                if row_norm <= scaler.config.zero_tolerance
                else 1.0 / row_norm
            )
        factor = min(
            scaler.config.max_constraint_scaling_factor,
            max(scaler.config.min_constraint_scaling_factor, factor),
        )
        set_scaling_factor(constraint, factor, overwrite=True)
    return {
        **cleared,
        "active_variables_scaled": len(active_variables),
        "active_constraints_scaled": len(active_constraints),
        "inherited_inequality_factors_preserved": len(
            inherited_inequality_factors
        ),
    }


def _checkpoint(name: str) -> Path:
    return REPO_ROOT / "initialization_files" / name


def _design_variables(model: Any) -> dict[str, Any]:
    """Return the eight design variables opened by the published optimization."""
    return {
        "Qs": model.fs.Qs,
        "H103_temperature": model.fs.H103.outlet.temperature[0],
        "H104_temperature": model.fs.H104.outlet.temperature[0],
        "H105_temperature": model.fs.H105.outlet.temperature[0],
        "F101_deltaP": model.fs.F101.deltaP[0],
        "H106_temperature": model.fs.H106.outlet.temperature[0],
        "H106_pressure": model.fs.H106.outlet.pressure[0],
        "F102_deltaP": model.fs.F102.deltaP[0],
    }


def _build_preoptimization_model(
    model_code: int = 5,
    region: str = "Bakken",
    costing_tax: float = 0.0,
    unit_initialization_region: str | None = None,
) -> Any:
    """Rebuild and load the published checkpoint chain before optimization."""
    model = create_flowsheet(model_code)
    define_models(model, catalyst_mass=1167.003367)
    define_arcs(model)

    inlet_data = pd.read_csv(REPO_ROOT / "data" / "NGL_compositions.csv")
    composition = {
        row["Species"]: (1e-6 if row[region] == 0.0 else round(row[region], 4))
        for _, row in inlet_data.iterrows()
    }
    conversion = {"ethane": 0.3566, "propane": 0.6632, "nbutane": 0.5188}
    set_unit_model_variables(
        model,
        model_code=model_code,
        feed_flow_rate=481.3888889,
        feed_temp=308.0,
        feed_pressure=700000.0,
        inlet_composition_dict=composition,
        dehydro_conv_dict=conversion,
    )
    set_scaling_factors(
        model,
        flow_mol_scaling_factor=1e-2 if model_code in (2, 3) else 1e-3,
        inlet_composition_dict=composition,
    )

    unit_region = unit_initialization_region or region
    ms.from_json(
        model,
        fname=str(
            _checkpoint(
                f"CISTAR_unit_initialization_{unit_region}_M{model_code}.json.gz"
            )
        ),
    )
    update_model_after_initialization(model)
    # Preserve the published notebook sequence exactly. The second call is
    # idempotent and reports that no degenerate constraints remain.
    vapor_only_to_vapor_liquid_reformulate(model.fs.T102)
    vapor_only_to_vapor_liquid_reformulate(model.fs.T102)
    ms.from_json(
        model,
        fname=str(
            _checkpoint(
                f"CISTAR_solve_constrained_{region}_M{model_code}_purge_0.01.json.gz"
            )
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
        region,
        str(REPO_ROOT / "data" / "LHV.xlsx"),
        str(REPO_ROOT / "data" / "NGL_compositions.csv"),
        str(REPO_ROOT / "data" / "NGL_fraction.csv"),
    )
    calculate_stream_energies(model)
    calculate_emissions(
        model,
        region,
        str(REPO_ROOT / "data" / "emissions_factor_by_region.csv"),
    )
    create_ghg_objective(model)
    calculate_costs_for_objective(model, c_tax_flag=True, c_tax_val=costing_tax)
    ms.from_json(
        model,
        fname=str(
            _checkpoint(
                f"CISTAR_solve_with_costing_{region}_C_tax_{costing_tax}_"
                f"M{model_code}_purge_0.01.json.gz"
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


def _write_report(report: dict[str, Any], output: Path | None) -> None:
    """Persist an incremental report without leaving a partially written JSON file."""
    if output is None:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(output)


def _load_column_order(model: Any, path: Path) -> tuple[list[Any], str]:
    """Resolve a recorded NL variable sequence against a rebuilt model."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    try:
        names = payload["ordering"]["variables"]["names"]
        expected_digest = payload["ordering"]["variables"]["sha256"]
    except (KeyError, TypeError) as err:
        raise ValueError(f"Not an NL symbol-map record: {path}") from err
    if not isinstance(names, list) or not all(isinstance(name, str) for name in names):
        raise ValueError(f"Variable ordering must be a list of names: {path}")
    digest = hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()
    if digest != expected_digest:
        raise ValueError(f"Variable-order digest does not match its names: {path}")
    components = []
    missing = []
    for name in names:
        component = model.find_component(name)
        if component is None:
            missing.append(name)
        else:
            components.append(component)
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(
            f"Column-order map has {len(missing)} unresolved variables: {preview}"
        )
    return components, digest


def _capture_named_state(
    model: Any,
    output: Path,
    *,
    stage: str,
    tax_rate: float,
) -> dict[str, Any]:
    """Write a name-aligned snapshot of every model variable."""
    variables = []
    missing_values = 0
    for variable in model.component_data_objects(
        Var, active=None, descend_into=True, sort=True
    ):
        variable_value = value(variable, exception=False)
        if variable_value is None:
            missing_values += 1
        variables.append(
            {
                "name": variable.name,
                "value": variable_value,
                "fixed": bool(variable.fixed),
                "lower_bound": value(variable.lb, exception=False),
                "upper_bound": value(variable.ub, exception=False),
            }
        )
    payload = {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-NAMED-STATE-SNAPSHOT-001",
        "stage": stage,
        "co2_tax_usd_per_kg": tax_rate,
        "environment": {
            "idaes": idaes.__version__,
            "pyomo": pyomo.__version__,
            "python": sys.version,
        },
        "variable_count": len(variables),
        "missing_value_count": missing_values,
        "variables": variables,
    }
    _write_report(payload, output)
    contents = output.read_bytes()
    return {
        "path": str(output.resolve()),
        "bytes": len(contents),
        "sha256": hashlib.sha256(contents).hexdigest(),
        "variable_count": len(variables),
        "missing_value_count": missing_values,
    }


def _tax_path_token(tax_rate: float) -> str:
    return format(tax_rate, ".12g").replace("-", "m").replace("+", "p")


def _regularize_pseudo_zero_inlet_phases(
    model: Any,
    threshold: float = 1e-6,
) -> dict[str, Any]:
    """Eliminate ill-conditioned composition equations for structural zero phases."""
    targets = (
        ("H105.inlet", model.fs.H105.control_volume.properties_in[0.0], "Liq"),
        ("H106.inlet", model.fs.H106.control_volume.properties_in[0.0], "Vap"),
        ("T103.inlet", model.fs.T103.properties_in[0.0], "Liq"),
        ("T104.inlet", model.fs.T104.properties_in[0.0], "Liq"),
    )
    details = []
    for label, state, phase in targets:
        components = [
            component
            for candidate_phase, component in state.phase_component_set
            if candidate_phase == phase
        ]
        flows = [
            max(0.0, float(value(state.flow_mol_phase_comp[phase, component])))
            for component in components
        ]
        total = sum(flows)
        if total >= threshold:
            raise RuntimeError(
                f"{label} {phase} flow {total} is not below {threshold}; "
                "refusing pseudo-zero-phase regularization."
            )
        fractions = (
            [flow / total for flow in flows]
            if total > 0.0
            else [1.0 / len(components)] * len(components)
        )
        deactivated = 0
        for component, fraction in zip(components, fractions):
            constraint = state.mole_frac_phase_comp_eq[phase, component]
            if constraint.active:
                constraint.deactivate()
                deactivated += 1
            state.mole_frac_phase_comp[phase, component].fix(fraction)
        details.append(
            {
                "state": label,
                "phase": phase,
                "component_count": len(components),
                "starting_total_phase_flow_mol_per_s": total,
                "constraints_deactivated": deactivated,
                "compositions_fixed": len(components),
                "fixed_fraction_sum": sum(fractions),
            }
        )
    return {
        "threshold_mol_per_s": threshold,
        "states_regularized": len(details),
        "constraints_deactivated": sum(
            item["constraints_deactivated"] for item in details
        ),
        "compositions_fixed": sum(item["compositions_fixed"] for item in details),
        "details": details,
    }


def _apply_r102_targeted_scaling(
    model: Any,
    heat_scaling_factor: float | None,
    normalize_rows: bool,
) -> dict[str, Any]:
    """Apply independently selectable R102 heat and row-norm scaling controls."""
    heat_variables = list(model.fs.R102.control_volume.heat.values())
    if heat_scaling_factor is not None:
        for variable in heat_variables:
            set_scaling_factor(variable, heat_scaling_factor, overwrite=True)

    report = {
        "heat_variables_scaled": (
            len(heat_variables) if heat_scaling_factor is not None else 0
        ),
        "heat_scaling_factor": heat_scaling_factor,
        "row_normalization_enabled": normalize_rows,
    }
    if not normalize_rows:
        report["constraints_row_normalized"] = 0
        return report

    jacobian, nlp = get_jacobian(model, scaled=True)
    row_norms = np.asarray(sparse.linalg.norm(jacobian, ord=2, axis=1)).ravel()
    selected_norms = []
    zero_norm_constraints = []
    for index, constraint in enumerate(nlp.clist):
        if not constraint.name.startswith("fs.R102"):
            continue
        norm = float(row_norms[index])
        if norm == 0.0:
            zero_norm_constraints.append(constraint.name)
            continue
        existing = get_constraint_transform_applied_scaling_factor(
            constraint, default=1.0
        )
        constraint_scaling_transform(
            constraint, float(existing) / norm, overwrite=True
        )
        selected_norms.append(norm)
    report.update({
        "constraints_row_normalized": len(selected_norms),
        "zero_norm_constraints_skipped": zero_norm_constraints,
        "pre_normalization_scaled_row_norm_2_min": min(selected_norms),
        "pre_normalization_scaled_row_norm_2_median": float(
            np.median(selected_norms)
        ),
        "pre_normalization_scaled_row_norm_2_max": max(selected_norms),
    })
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ipopt", type=Path, required=True)
    parser.add_argument(
        "--linear-solver", choices=("ma27", "ma57"), default="ma27"
    )
    parser.add_argument(
        "--ma57-automatic-scaling",
        action="store_true",
        help=(
            "Enable MA57's internal automatic scaling (ICNTL(15)); valid only "
            "with --linear-solver ma57."
        ),
    )
    parser.add_argument(
        "--regularize-pseudo-zero-inlet-phases",
        action="store_true",
        help=(
            "Fix analytically normalized compositions and deactivate their "
            "ill-conditioned normalization equations for four structurally "
            "absent inlet phases (H105 Liq, H106 Vap, T103 Liq, T104 Liq)."
        ),
    )
    parser.add_argument(
        "--r102-heat-scaling-factor",
        type=float,
        help=(
            "Override scaling factors for R102 heat variables. A magnitude-"
            "based diagnostic value is 1e-8. Requires --nlp-scaling-method "
            "user-scaling."
        ),
    )
    parser.add_argument(
        "--r102-row-norm-scaling",
        action="store_true",
        help=(
            "Normalize R102 constraint rows by their current scaled-Jacobian "
            "2-norm. Requires --nlp-scaling-method user-scaling."
        ),
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
        "--acceptable-tol",
        type=float,
        help=(
            "Opt in to IPOPT's acceptable-termination criterion at this KKT "
            "tolerance. The strict tol remains 1e-6. Intended only for "
            "diagnosing paths that cross an accurate iterate before "
            "restoration failure."
        ),
    )
    parser.add_argument(
        "--acceptable-iter",
        type=int,
        default=1,
        help="Consecutive acceptable iterates required when --acceptable-tol is set.",
    )
    parser.add_argument(
        "--nlp-scaling-method",
        choices=("gradient-based", "user-scaling", "none"),
        help="Set Ipopt's nlp_scaling_method; omit to retain its default.",
    )
    scaling_group = parser.add_mutually_exclusive_group()
    scaling_group.add_argument(
        "--modern-autoscale",
        action="store_true",
        help=(
            "Overwrite variable scaling by current magnitude and constraint "
            "scaling by Jacobian norm using the IDAES 2.12 AutoScaler."
        ),
    )
    scaling_group.add_argument(
        "--active-nlp-autoscale",
        action="store_true",
        help=(
            "Clear inherited scaling suffixes, scale active NLP variables by "
            "magnitude, and scale active constraints by Jacobian norm."
        ),
    )
    parser.add_argument("--autoscale-norm", type=int, default=2)
    parser.add_argument(
        "--preserve-inherited-inequality-scaling",
        action="store_true",
        help=(
            "With --active-nlp-autoscale, retain inherited factors for active "
            "inequalities while rebuilding variable and equality factors."
        ),
    )
    parser.add_argument(
        "--solver-io",
        choices=("nl", "nl_v1", "nl_v2"),
        default="nl",
        help="Select the Pyomo AMPL NL writer implementation.",
    )
    parser.add_argument(
        "--inline-defined-variables",
        action="store_true",
        help=(
            "Set the modern NL writer's export_defined_variables option to "
            "false, inlining named Expression objects instead of exporting "
            "them as NL common expressions. Valid only with --solver-io nl."
        ),
    )
    parser.add_argument(
        "--file-determinism",
        choices=("ordered", "sort-indices", "sort-symbols"),
        default="ordered",
        help=(
            "Control deterministic NL row/column ordering. sort-symbols gives "
            "both environments a component-name-based ordering rule."
        ),
    )
    parser.add_argument(
        "--column-order-from-symbol-map",
        type=Path,
        help=(
            "Resolve and pass the variable sequence from an exported NL "
            "symbol-map JSON as Pyomo's explicit column_order. This is a "
            "cross-version diagnostic control, not a recommended default."
        ),
    )
    parser.add_argument(
        "--initial-optimal-tax",
        type=float,
        metavar="USD_PER_KG",
        help=(
            "Load the archived M5/Bakken optimum at this tax before the first "
            "requested solve. Useful for isolating one sequential transition."
        ),
    )
    parser.add_argument(
        "--free-design-variables",
        nargs="+",
        choices=DESIGN_VARIABLE_NAMES,
        help=(
            "Free only this subset of the eight published optimization "
            "variables; all omitted design variables remain at the checkpoint."
        ),
    )
    parser.add_argument("--tee", action="store_true")
    parser.add_argument(
        "--continue-after-failure",
        action="store_true",
        help=(
            "Continue a tax sequence after a non-optimal solver termination. "
            "By default the runner stops so a failed state is not used to "
            "initialize later continuation points."
        ),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--state-snapshot-dir",
        type=Path,
        help=(
            "Write complete named-variable snapshots immediately before and "
            "after each solve. Intended for detailed state-drift diagnostics."
        ),
    )
    args = parser.parse_args()
    if (
        args.preserve_inherited_inequality_scaling
        and not args.active_nlp_autoscale
    ):
        parser.error(
            "--preserve-inherited-inequality-scaling requires "
            "--active-nlp-autoscale"
        )
    if args.ma57_automatic_scaling and args.linear_solver != "ma57":
        parser.error("--ma57-automatic-scaling requires --linear-solver ma57")
    if args.inline_defined_variables and args.solver_io != "nl":
        parser.error("--inline-defined-variables requires --solver-io nl")
    if (
        args.r102_heat_scaling_factor is not None or args.r102_row_norm_scaling
    ) and args.nlp_scaling_method != "user-scaling":
        parser.error(
            "R102 scaling controls require --nlp-scaling-method user-scaling"
        )

    started = time.time()
    ipopt = args.ipopt.resolve()
    solver_environment = _configure_solver_environment(ipopt)
    report: dict[str, Any] = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            **_idaes_provenance(),
            "pyomo": pyomo.__version__,
            "ipopt": _solver_version(ipopt, solver_environment),
            "linear_solver": args.linear_solver,
            "ma57_automatic_scaling": args.ma57_automatic_scaling,
            "acceptable_tol": args.acceptable_tol,
            "acceptable_iter": (
                args.acceptable_iter if args.acceptable_tol is not None else None
            ),
            "acceptable_constr_viol_tol": (
                1e-6 if args.acceptable_tol is not None else None
            ),
            "acceptable_dual_inf_tol": args.acceptable_tol,
            "acceptable_compl_inf_tol": args.acceptable_tol,
            "regularize_pseudo_zero_inlet_phases": (
                args.regularize_pseudo_zero_inlet_phases
            ),
            "r102_heat_scaling_factor": args.r102_heat_scaling_factor,
            "r102_row_norm_scaling": args.r102_row_norm_scaling,
            "nlp_scaling_method": args.nlp_scaling_method or "ipopt-default",
            "modern_autoscale": args.modern_autoscale,
            "active_nlp_autoscale": args.active_nlp_autoscale,
            "preserve_inherited_inequality_scaling": (
                args.preserve_inherited_inequality_scaling
            ),
            "autoscale_norm": (
                args.autoscale_norm
                if args.modern_autoscale or args.active_nlp_autoscale
                else None
            ),
            "solver_io": args.solver_io,
            "inline_defined_variables": args.inline_defined_variables,
            "continue_after_failure": args.continue_after_failure,
            "file_determinism": args.file_determinism,
            "column_order_from_symbol_map": (
                str(args.column_order_from_symbol_map.resolve())
                if args.column_order_from_symbol_map is not None
                else None
            ),
        },
        "case": {"model_code": 5, "region": "Bakken"},
        "free_design_variables": (
            list(args.free_design_variables)
            if args.free_design_variables is not None
            else list(DESIGN_VARIABLE_NAMES)
        ),
        "tax_rate_units": "USD/kg CO2e",
        "initialization_chain": [
            "CISTAR_unit_initialization_Bakken_M5.json.gz",
            "CISTAR_solve_constrained_Bakken_M5_purge_0.01.json.gz",
            "CISTAR_solve_with_costing_Bakken_C_tax_0.0_M5_purge_0.01.json.gz",
        ],
        "status": "building",
        "runs": [],
    }

    model = _build_preoptimization_model()
    column_order = None
    if args.column_order_from_symbol_map is not None:
        column_order, column_order_digest = _load_column_order(
            model, args.column_order_from_symbol_map
        )
        report["environment"]["requested_column_order_sha256"] = (
            column_order_digest
        )
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
    pseudo_zero_phase_regularization = None
    if args.regularize_pseudo_zero_inlet_phases:
        pseudo_zero_phase_regularization = _regularize_pseudo_zero_inlet_phases(
            model
        )
        report["pseudo_zero_phase_regularization"] = (
            pseudo_zero_phase_regularization
        )
    report["build_seconds"] = time.time() - started
    report["status"] = "running"
    report["total_wall_seconds"] = time.time() - started
    _write_report(report, args.output)
    archived_rows = _load_archived_rows()
    r102_targeted_scaling = None

    for tax_rate in args.tax_rates:
        case_started = time.time()
        model.fs.c_tax_rate = tax_rate
        unfix_DOFs_pre_optimization(model)
        if (
            args.r102_heat_scaling_factor is not None
            or args.r102_row_norm_scaling
        ) and r102_targeted_scaling is None:
            r102_targeted_scaling = _apply_r102_targeted_scaling(
                model,
                args.r102_heat_scaling_factor,
                args.r102_row_norm_scaling,
            )
            report["r102_targeted_scaling"] = r102_targeted_scaling
        if args.free_design_variables is not None:
            selected = set(args.free_design_variables)
            for name, variable in _design_variables(model).items():
                if name not in selected:
                    variable.fix()
        initial_dof = degrees_of_freedom(model)
        starting_results = _collect_results(model)
        state_snapshots = None
        if args.state_snapshot_dir is not None:
            token = _tax_path_token(tax_rate)
            state_snapshots = {
                "start": _capture_named_state(
                    model,
                    args.state_snapshot_dir / f"tax-{token}-start.json",
                    stage="solver-start",
                    tax_rate=tax_rate,
                )
            }
        archived = archived_rows.get(tax_rate * 1000)
        starting_comparison = None
        if archived is not None:
            starting_comparison = {
                metric: starting_results[metric] - archived[metric]
                for metric in starting_results
            }
        if args.modern_autoscale:
            _apply_modern_autoscaling(model, args.autoscale_norm)
        active_nlp_scaling = None
        if args.active_nlp_autoscale:
            active_nlp_scaling = _apply_active_nlp_autoscaling(
                model,
                args.autoscale_norm,
                preserve_inherited_inequality_scaling=(
                    args.preserve_inherited_inequality_scaling
                ),
            )
        solver = SolverFactory(
            "ipopt", solver_io=args.solver_io, executable=str(ipopt)
        )
        solver.options.update(
            {
                "tol": 1e-6,
                "bound_push": 1e-8,
                "max_iter": args.max_iter,
                "linear_solver": args.linear_solver,
            }
        )
        if args.nlp_scaling_method is not None:
            solver.options["nlp_scaling_method"] = args.nlp_scaling_method
        if args.ma57_automatic_scaling:
            solver.options["ma57_automatic_scaling"] = "yes"
        if args.acceptable_tol is not None:
            solver.options["acceptable_tol"] = args.acceptable_tol
            solver.options["acceptable_iter"] = args.acceptable_iter
            # IPOPT's aggregate acceptable_tol uses scaled quantities, while
            # the component thresholds otherwise have very loose defaults.
            # Bind them explicitly so an "acceptable" diagnostic cannot hide
            # a large unscaled residual or dual infeasibility.
            solver.options["acceptable_constr_viol_tol"] = 1e-6
            solver.options["acceptable_dual_inf_tol"] = args.acceptable_tol
            solver.options["acceptable_compl_inf_tol"] = args.acceptable_tol
        determinism = {
            "ordered": FileDeterminism.ORDERED,
            "sort-indices": FileDeterminism.SORT_INDICES,
            "sort-symbols": FileDeterminism.SORT_SYMBOLS,
        }[args.file_determinism]
        writer_options = {"file_determinism": determinism}
        if args.inline_defined_variables:
            writer_options["export_defined_variables"] = False
        if column_order is not None:
            writer_options["column_order"] = column_order
        solve_result = solver.solve(
            model, tee=args.tee, load_solutions=False, **writer_options
        )
        load_error = None
        try:
            model.solutions.load_from(solve_result)
        except ValueError as err:
            load_error = str(err)
        fresh = _collect_results(model)
        if state_snapshots is not None:
            state_snapshots["final"] = _capture_named_state(
                model,
                args.state_snapshot_dir / f"tax-{token}-final.json",
                stage="solver-final",
                tax_rate=tax_rate,
            )
        comparison = None
        if archived is not None:
            comparison = {
                metric: fresh[metric] - archived[metric] for metric in fresh
            }
        report["runs"].append(
            {
                "co2_tax_usd_per_kg": tax_rate,
                "initial_degrees_of_freedom": initial_dof,
                "starting_results_in_migrated_csv_units": starting_results,
                "starting_difference_from_migrated_csv": starting_comparison,
                "termination_condition": str(
                    solve_result.solver.termination_condition
                ),
                "solver_status": str(solve_result.solver.status),
                "solver_message": str(solve_result.solver.message),
                "solution_load_error": load_error,
                "active_nlp_scaling": active_nlp_scaling,
                "pseudo_zero_phase_regularization": (
                    pseudo_zero_phase_regularization
                ),
                "r102_targeted_scaling": r102_targeted_scaling,
                "state_snapshots": state_snapshots,
                "wall_seconds": time.time() - case_started,
                "results_in_migrated_csv_units": fresh,
                "difference_from_migrated_csv": comparison,
                "figure_data": collect_figure_data(model),
            }
        )
        fix_DOFs_post_optimization(model)
        report["total_wall_seconds"] = time.time() - started
        _write_report(report, args.output)
        if (
            str(solve_result.solver.termination_condition) != "optimal"
            and not args.continue_after_failure
        ):
            report["terminated_early_after_solver_failure"] = True
            break

    report["total_wall_seconds"] = time.time() - started
    report["status"] = "complete"
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0 if all(
        run["termination_condition"] == "optimal" for run in report["runs"]
    ) else 2


if __name__ == "__main__":
    raise SystemExit(main())
