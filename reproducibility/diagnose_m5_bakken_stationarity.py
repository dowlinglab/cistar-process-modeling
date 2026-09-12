#!/usr/bin/env python3
"""Measure first-order stationarity at the archived M5/Bakken optimum."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsmr, splu

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.scaling import get_jacobian, get_scaling_factor
from pyomo.environ import value

from reproducibility.run_m5_bakken_tax_series import (
    _apply_modern_autoscaling,
    _build_preoptimization_model,
    _checkpoint,
    _collect_results,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


ARCHIVED_OPTIMUM = (
    "CISTAR_optimal_solution_Bakken_C_tax_0.0_M5_purge_0.01_"
    "sequential_solve.json.gz"
)


def stationarity_projection(
    gradient: np.ndarray,
    gradient_rows: sparse.spmatrix,
    *,
    atol: float = 1e-10,
    btol: float = 1e-10,
    maxiter: int = 20000,
) -> dict[str, Any]:
    """Project a gradient onto a sparse row span and summarize the residual."""
    matrix = sparse.csr_matrix(gradient_rows, dtype=float)
    gradient = np.asarray(gradient, dtype=float)
    solution = lsmr(
        matrix.transpose(),
        -gradient,
        atol=atol,
        btol=btol,
        maxiter=maxiter,
    )
    residual = gradient + matrix.transpose().dot(solution[0])
    gradient_norm = float(np.linalg.norm(gradient))
    residual_norm = float(np.linalg.norm(residual))
    return {
        "gradient_norm_2": gradient_norm,
        "stationarity_residual_norm_2": residual_norm,
        "relative_stationarity_residual": (
            residual_norm / gradient_norm if gradient_norm else 0.0
        ),
        "lsmr_stop_code": int(solution[1]),
        "lsmr_iterations": int(solution[2]),
        "estimated_condition": float(solution[6]),
        "residual": residual,
    }


def _active_rows(
    jacobian: sparse.csr_matrix,
    constraints: list[Any],
    tolerance: float,
) -> tuple[sparse.csr_matrix, list[str], list[str], list[int]]:
    indices: list[int] = []
    names: list[str] = []
    inequality_names: list[str] = []
    inequality_indices: list[int] = []
    for index, constraint in enumerate(constraints):
        if constraint.equality:
            indices.append(index)
            names.append(constraint.name)
            continue
        body = float(value(constraint.body))
        lower = (
            float(value(constraint.lower))
            if constraint.has_lb()
            else None
        )
        upper = (
            float(value(constraint.upper))
            if constraint.has_ub()
            else None
        )
        active = (
            lower is not None
            and abs(body - lower) <= tolerance * max(1.0, abs(lower))
        ) or (
            upper is not None
            and abs(body - upper) <= tolerance * max(1.0, abs(upper))
        )
        if active:
            indices.append(index)
            names.append(constraint.name)
            inequality_names.append(constraint.name)
            inequality_indices.append(index)
    return jacobian[indices, :], names, inequality_names, inequality_indices


def reduced_space_gradient(
    gradient: np.ndarray,
    equality_jacobian: sparse.spmatrix,
    design_indices: list[int],
    active_inequality_jacobian: sparse.spmatrix | None = None,
) -> dict[str, Any]:
    """Eliminate square state equations and return design-space gradients."""
    jacobian = sparse.csr_matrix(equality_jacobian, dtype=float)
    design_set = set(design_indices)
    state_indices = [
        index for index in range(jacobian.shape[1]) if index not in design_set
    ]
    state_jacobian = jacobian[:, state_indices].tocsc()
    if state_jacobian.shape[0] != state_jacobian.shape[1]:
        raise ValueError(
            "State Jacobian is not square after removing the design variables: "
            f"{state_jacobian.shape}."
        )
    factorization = splu(state_jacobian)
    design_jacobian = jacobian[:, design_indices].toarray()
    state_sensitivities = factorization.solve(-design_jacobian)
    gradient = np.asarray(gradient, dtype=float)
    reduced_objective = (
        gradient[design_indices]
        + state_sensitivities.transpose().dot(gradient[state_indices])
    )
    result: dict[str, Any] = {
        "state_variables": len(state_indices),
        "design_variables": len(design_indices),
        "reduced_objective_gradient": reduced_objective,
        "reduced_objective_gradient_norm_2": float(
            np.linalg.norm(reduced_objective)
        ),
        "state_jacobian_lu_u_diagonal_ratio": float(
            np.max(np.abs(factorization.U.diagonal()))
            / np.min(np.abs(factorization.U.diagonal()))
        ),
    }
    if active_inequality_jacobian is not None:
        inequality = sparse.csr_matrix(active_inequality_jacobian, dtype=float)
        reduced_inequality = (
            inequality[:, design_indices].toarray()
            + inequality[:, state_indices].dot(state_sensitivities)
        )
        result["reduced_active_inequality_jacobian"] = reduced_inequality
    return result


def _active_bound_rows(
    variables: list[Any], tolerance: float
) -> tuple[sparse.csr_matrix, list[str]]:
    indices: list[int] = []
    names: list[str] = []
    for index, variable in enumerate(variables):
        current = float(value(variable))
        lower = float(value(variable.lb)) if variable.lb is not None else None
        upper = float(value(variable.ub)) if variable.ub is not None else None
        active_lower = lower is not None and abs(current - lower) <= (
            tolerance * max(1.0, abs(lower))
        )
        active_upper = upper is not None and abs(current - upper) <= (
            tolerance * max(1.0, abs(upper))
        )
        if active_lower or active_upper:
            indices.append(index)
            side = "lower" if active_lower else "upper"
            names.append(f"{variable.name}:{side}")
    rows = sparse.csr_matrix(
        (
            np.ones(len(indices)),
            (np.arange(len(indices), dtype=int), np.asarray(indices, dtype=int)),
        ),
        shape=(len(indices), len(variables)),
    )
    return rows, names


def diagnose(
    scaling: str,
    active_tolerance: float,
    top: int,
) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    model.fs.c_tax_rate = 0.0
    unfix_DOFs_pre_optimization(model)
    checkpoint = _checkpoint(ARCHIVED_OPTIMUM)
    ms.from_json(model, fname=str(checkpoint))
    # The serialized optimum restores the historical fixed flags, so reopen the
    # eight published design degrees of freedom after loading it.
    unfix_DOFs_pre_optimization(model)

    if scaling == "modern-auto":
        _apply_modern_autoscaling(model, norm=2)
    elif scaling != "legacy":
        raise ValueError(f"Unsupported scaling mode: {scaling}")

    jacobian, nlp = get_jacobian(model, scaled=True)
    variables = nlp.vlist
    constraints = nlp.clist
    variable_scaling = np.asarray(
        [get_scaling_factor(variable, default=1.0) for variable in variables]
    )
    scaled_gradient = nlp.evaluate_grad_objective() / variable_scaling

    equality_indices = [
        index for index, constraint in enumerate(constraints) if constraint.equality
    ]
    equality_rows = jacobian[equality_indices, :]
    (
        active_constraint_rows,
        active_constraint_names,
        active_inequality_names,
        active_inequality_indices,
    ) = _active_rows(jacobian, constraints, active_tolerance)
    bound_rows, active_bound_names = _active_bound_rows(
        variables, active_tolerance
    )
    equality_projection = stationarity_projection(scaled_gradient, equality_rows)
    active_projection = stationarity_projection(
        scaled_gradient,
        sparse.vstack([active_constraint_rows, bound_rows], format="csr"),
    )

    residual = active_projection.pop("residual")
    equality_projection.pop("residual")
    ranked = np.argsort(np.abs(residual))[::-1][:top]

    design_objects = [
        model.fs.Qs,
        model.fs.H103.outlet.temperature[0],
        model.fs.H104.outlet.temperature[0],
        model.fs.H105.outlet.temperature[0],
        model.fs.F101.deltaP[0],
        model.fs.H106.outlet.temperature[0],
        model.fs.H106.outlet.pressure[0],
        model.fs.F102.deltaP[0],
    ]
    variable_index = {variable.name: index for index, variable in enumerate(variables)}
    design_indices = [variable_index[variable.name] for variable in design_objects]
    reduced = reduced_space_gradient(
        scaled_gradient,
        equality_rows,
        design_indices,
        jacobian[active_inequality_indices, :],
    )
    reduced_objective = reduced.pop("reduced_objective_gradient")
    reduced_inequality = reduced.pop("reduced_active_inequality_jacobian")
    reduced["objective_gradient_by_design_variable"] = [
        {
            "variable": variable.name,
            "value": float(value(variable)),
            "scaling_factor": float(variable_scaling[index]),
            "reduced_gradient": float(component),
        }
        for variable, index, component in zip(
            design_objects, design_indices, reduced_objective
        )
    ]
    reduced["active_inequality_reduced_gradients"] = [
        {
            "constraint": name,
            "gradient": [float(component) for component in row],
        }
        for name, row in zip(active_inequality_names, reduced_inequality)
    ]
    if len(active_inequality_names):
        multiplier, _, _, _ = np.linalg.lstsq(
            reduced_inequality.transpose(), -reduced_objective, rcond=None
        )
        reduced_residual = (
            reduced_objective + reduced_inequality.transpose().dot(multiplier)
        )
    else:
        multiplier = np.empty(0)
        reduced_residual = reduced_objective
    reduced["unconstrained_active_inequality_multipliers"] = [
        {"constraint": name, "multiplier": float(component)}
        for name, component in zip(active_inequality_names, multiplier)
    ]
    reduced["active_set_stationarity_residual"] = [
        float(component) for component in reduced_residual
    ]
    reduced["active_set_stationarity_residual_norm_2"] = float(
        np.linalg.norm(reduced_residual)
    )
    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-ARCHIVED-OPTIMUM-KKT-001",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "pyomo": pyomo.__version__,
        },
        "case": {
            "model_code": 5,
            "region": "Bakken",
            "co2_tax_usd_per_kg": 0.0,
            "checkpoint": checkpoint.name,
            "scaling": scaling,
            "active_tolerance": active_tolerance,
        },
        "checkpoint_results_in_migrated_csv_units": _collect_results(model),
        "nlp": {
            "variables": len(variables),
            "constraints": len(constraints),
            "equalities": len(equality_indices),
            "active_constraints_including_equalities": len(active_constraint_names),
            "active_inequalities": len(active_inequality_names),
            "active_variable_bounds": len(active_bound_names),
            "jacobian_nonzeros": int(jacobian.nnz),
        },
        "equality_only_projection": equality_projection,
        "active_set_projection": active_projection,
        "reduced_space_analysis": reduced,
        "active_inequality_names": active_inequality_names,
        "active_bound_names": active_bound_names,
        "largest_unexplained_scaled_gradient_components": [
            {
                "variable": variables[index].name,
                "residual": float(residual[index]),
                "absolute_residual": float(abs(residual[index])),
                "value": float(value(variables[index])),
                "scaling_factor": float(variable_scaling[index]),
            }
            for index in ranked
        ],
        "interpretation_limit": (
            "The active-set projection allows multipliers of either sign. A small "
            "residual is necessary but not sufficient for full KKT optimality; a "
            "large residual proves that the archived point is not stationary for "
            "the represented active-set gradients at the requested tolerance."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaling", choices=("legacy", "modern-auto"), default="legacy")
    parser.add_argument("--active-tolerance", type=float, default=1e-7)
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = diagnose(args.scaling, args.active_tolerance, args.top)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
