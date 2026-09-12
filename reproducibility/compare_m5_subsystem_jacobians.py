#!/usr/bin/env python3
"""Compare named subsystem Jacobian rows and columns at two M5 states."""

from __future__ import annotations

import argparse
import json
import math
import platform
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.scaling import get_jacobian, get_scaling_factor
from pyomo.environ import Var, value

from reproducibility.run_m5_bakken_tax_series import (
    _build_preoptimization_model,
    _checkpoint,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


ARCHIVED_OPTIMUM = (
    "CISTAR_optimal_solution_Bakken_C_tax_0.0_M5_purge_0.01_"
    "sequential_solve.json.gz"
)


def _axis_metrics(matrix: sparse.spmatrix, axis: int) -> dict[str, np.ndarray]:
    csr = sparse.csr_matrix(matrix, dtype=float)
    norms = np.asarray(sparse.linalg.norm(csr, ord=2, axis=axis)).ravel()
    maxima = np.asarray(abs(csr).max(axis=axis).toarray()).ravel()
    nonzeros = (
        np.diff(csr.indptr)
        if axis == 1
        else np.diff(csr.tocsc().indptr)
    )
    return {"norm_2": norms, "max_abs": maxima, "nonzeros": nonzeros}


def _change_ratio(reference: float, candidate: float) -> float | None:
    if reference == candidate == 0.0:
        return 1.0
    smaller = min(abs(reference), abs(candidate))
    if smaller == 0.0:
        return None
    return max(abs(reference), abs(candidate)) / smaller


def _component_family(name: str) -> str:
    return re.sub(r"\[[^]]*\]", "", name)


def _family_summaries(
    records: list[dict[str, Any]], norm_key: str
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[_component_family(record["name"])].append(record)
    summaries = []
    for family, members in grouped.items():
        norms = np.asarray([member[norm_key] for member in members], dtype=float)
        summaries.append(
            {
                "family": family,
                "count": len(members),
                "scaled_norm_2_min": float(np.min(norms)),
                "scaled_norm_2_median": float(np.median(norms)),
                "scaled_norm_2_max": float(np.max(norms)),
                "scaling_factors": sorted(
                    {member["scaling_factor"] for member in members}
                ),
            }
        )
    return sorted(
        summaries, key=lambda item: item["scaled_norm_2_max"], reverse=True
    )


def _apply_snapshot(model: Any, snapshot_path: Path) -> dict[str, int]:
    payload = json.loads(snapshot_path.read_text())
    variables = {
        variable.name: variable
        for variable in model.component_data_objects(
            Var, active=None, descend_into=True
        )
    }
    applied = 0
    missing_in_model = 0
    missing_values = 0
    for item in payload["variables"]:
        variable = variables.get(item["name"])
        if variable is None:
            missing_in_model += 1
            continue
        if item["value"] is None:
            missing_values += 1
            continue
        variable.set_value(item["value"], skip_validation=True)
        applied += 1
    return {
        "snapshot_variables": len(payload["variables"]),
        "values_applied": applied,
        "names_missing_in_model": missing_in_model,
        "missing_values_skipped": missing_values,
    }


def _capture(model: Any) -> dict[str, Any]:
    unscaled, unscaled_nlp = get_jacobian(model, scaled=False)
    scaled, scaled_nlp = get_jacobian(model, scaled=True)
    variable_names = [variable.name for variable in unscaled_nlp.vlist]
    constraint_names = [constraint.name for constraint in unscaled_nlp.clist]
    if variable_names != [variable.name for variable in scaled_nlp.vlist]:
        raise RuntimeError("Scaled and unscaled variable orders differ.")
    if constraint_names != [constraint.name for constraint in scaled_nlp.clist]:
        raise RuntimeError("Scaled and unscaled constraint orders differ.")
    return {
        "variable_names": variable_names,
        "constraint_names": constraint_names,
        "variables": unscaled_nlp.vlist,
        "constraints": unscaled_nlp.clist,
        "variable_values": [
            _finite_value(variable) for variable in unscaled_nlp.vlist
        ],
        "unscaled_columns": _axis_metrics(unscaled, axis=0),
        "scaled_columns": _axis_metrics(scaled, axis=0),
        "unscaled_rows": _axis_metrics(unscaled, axis=1),
        "scaled_rows": _axis_metrics(scaled, axis=1),
    }


def _finite_value(component: Any) -> float | None:
    result = value(component, exception=False)
    return float(result) if result is not None and math.isfinite(result) else None


def compare(snapshot_path: Path, prefix: str, top: int) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    model.fs.c_tax_rate = 0.0
    unfix_DOFs_pre_optimization(model)
    checkpoint = _checkpoint(ARCHIVED_OPTIMUM)
    ms.from_json(model, fname=str(checkpoint))
    unfix_DOFs_pre_optimization(model)
    reference = _capture(model)

    application = _apply_snapshot(model, snapshot_path)
    candidate = _capture(model)
    if reference["variable_names"] != candidate["variable_names"]:
        raise RuntimeError("Active variable names differ between states.")
    if reference["constraint_names"] != candidate["constraint_names"]:
        raise RuntimeError("Active constraint names differ between states.")

    variable_records = []
    for index, name in enumerate(reference["variable_names"]):
        if not name.startswith(prefix):
            continue
        variable = candidate["variables"][index]
        variable_records.append(
            {
                "name": name,
                "reference_value": reference["variable_values"][index],
                "candidate_value": _finite_value(variable),
                "scaling_factor": float(
                    get_scaling_factor(variable, default=1.0)
                ),
                "unscaled_column_norm_2_reference": float(
                    reference["unscaled_columns"]["norm_2"][index]
                ),
                "unscaled_column_norm_2_candidate": float(
                    candidate["unscaled_columns"]["norm_2"][index]
                ),
                "scaled_column_norm_2_reference": float(
                    reference["scaled_columns"]["norm_2"][index]
                ),
                "scaled_column_norm_2_candidate": float(
                    candidate["scaled_columns"]["norm_2"][index]
                ),
                "scaled_column_norm_change_ratio": _change_ratio(
                    reference["scaled_columns"]["norm_2"][index],
                    candidate["scaled_columns"]["norm_2"][index],
                ),
                "column_nonzeros": int(
                    candidate["scaled_columns"]["nonzeros"][index]
                ),
            }
        )
        if variable_records[-1]["candidate_value"] is not None:
            variable_records[-1]["candidate_scaled_absolute_value"] = abs(
                variable_records[-1]["candidate_value"]
                * variable_records[-1]["scaling_factor"]
            )

    constraint_records = []
    for index, name in enumerate(reference["constraint_names"]):
        if not name.startswith(prefix):
            continue
        constraint = candidate["constraints"][index]
        constraint_records.append(
            {
                "name": name,
                "scaling_factor": float(
                    get_scaling_factor(constraint, default=1.0)
                ),
                "unscaled_row_norm_2_reference": float(
                    reference["unscaled_rows"]["norm_2"][index]
                ),
                "unscaled_row_norm_2_candidate": float(
                    candidate["unscaled_rows"]["norm_2"][index]
                ),
                "scaled_row_norm_2_reference": float(
                    reference["scaled_rows"]["norm_2"][index]
                ),
                "scaled_row_norm_2_candidate": float(
                    candidate["scaled_rows"]["norm_2"][index]
                ),
                "scaled_row_norm_change_ratio": _change_ratio(
                    reference["scaled_rows"]["norm_2"][index],
                    candidate["scaled_rows"]["norm_2"][index],
                ),
                "row_nonzeros": int(candidate["scaled_rows"]["nonzeros"][index]),
            }
        )

    def largest(records: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
        return sorted(
            records,
            key=lambda item: (
                float("inf") if item[key] is None else item[key]
            ),
            reverse=True,
        )[:top]

    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-SUBSYSTEM-JACOBIAN-COMPARISON-001",
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
            "reference_checkpoint": checkpoint.name,
            "candidate_snapshot": str(snapshot_path),
            "component_prefix": prefix,
        },
        "snapshot_application": application,
        "active_nlp": {
            "variables": len(reference["variable_names"]),
            "constraints": len(reference["constraint_names"]),
            "subsystem_variables": len(variable_records),
            "subsystem_constraints": len(constraint_records),
        },
        "largest_scaled_column_norm_changes": largest(
            variable_records, "scaled_column_norm_change_ratio"
        ),
        "largest_scaled_row_norm_changes": largest(
            constraint_records, "scaled_row_norm_change_ratio"
        ),
        "largest_candidate_scaled_column_norms": sorted(
            variable_records,
            key=lambda item: item["scaled_column_norm_2_candidate"],
            reverse=True,
        )[:top],
        "smallest_nonzero_candidate_scaled_column_norms": sorted(
            (
                item
                for item in variable_records
                if item["scaled_column_norm_2_candidate"] > 0.0
            ),
            key=lambda item: item["scaled_column_norm_2_candidate"],
        )[:top],
        "largest_candidate_scaled_row_norms": sorted(
            constraint_records,
            key=lambda item: item["scaled_row_norm_2_candidate"],
            reverse=True,
        )[:top],
        "smallest_nonzero_candidate_scaled_row_norms": sorted(
            (
                item
                for item in constraint_records
                if item["scaled_row_norm_2_candidate"] > 0.0
            ),
            key=lambda item: item["scaled_row_norm_2_candidate"],
        )[:top],
        "variable_family_summaries": _family_summaries(
            variable_records, "scaled_column_norm_2_candidate"
        ),
        "constraint_family_summaries": _family_summaries(
            constraint_records, "scaled_row_norm_2_candidate"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-snapshot", type=Path, required=True)
    parser.add_argument("--prefix", default="fs.R102")
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = compare(arguments.candidate_snapshot, arguments.prefix, arguments.top)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    _write_report(report, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
