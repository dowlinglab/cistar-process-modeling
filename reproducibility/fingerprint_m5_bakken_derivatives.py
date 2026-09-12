#!/usr/bin/env python3
"""Fingerprint M5/Bakken derivatives at a shared archived or perturbed state."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from idaes.core.util import model_serializer as ms
from idaes.core.util.scaling import get_jacobian
from pyomo.environ import value

from reproducibility.diagnose_m5_bakken_stationarity import ARCHIVED_OPTIMUM
from reproducibility.run_m5_bakken_tax_series import (
    _build_preoptimization_model,
    _checkpoint,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


def canonical_sparse_fingerprints(
    matrix: sparse.spmatrix,
    row_names: list[str],
    column_names: list[str],
    significant_digits: Iterable[int] = (8, 10, 12),
) -> dict[str, Any]:
    """Hash sparse entries by component names, independent of storage order."""
    cleaned = sparse.csr_matrix(matrix, dtype=float)
    cleaned.sum_duplicates()
    cleaned.eliminate_zeros()
    coordinate = cleaned.tocoo()
    entries = sorted(
        (
            row_names[row],
            column_names[column],
            float(entry),
        )
        for row, column, entry in zip(
            coordinate.row, coordinate.col, coordinate.data
        )
    )

    def digest(render) -> str:
        hasher = hashlib.sha256()
        for row, column, entry in entries:
            hasher.update(row.encode("utf-8"))
            hasher.update(b"\0")
            hasher.update(column.encode("utf-8"))
            hasher.update(b"\0")
            hasher.update(render(entry).encode("ascii"))
            hasher.update(b"\n")
        return hasher.hexdigest()

    absolute_values = np.abs(coordinate.data)
    return {
        "shape": list(cleaned.shape),
        "nonzeros": int(cleaned.nnz),
        "frobenius_norm": float(np.linalg.norm(coordinate.data)),
        "minimum_nonzero_absolute_entry": float(np.min(absolute_values)),
        "maximum_absolute_entry": float(np.max(absolute_values)),
        "exact_float_hex_sha256": digest(float.hex),
        "rounded_sha256": {
            str(digits): digest(lambda entry, d=digits: format(entry, f".{d}g"))
            for digits in significant_digits
        },
    }


def canonical_vector_fingerprints(
    vector: np.ndarray,
    names: list[str],
    significant_digits: Iterable[int] = (8, 10, 12),
) -> dict[str, Any]:
    matrix = sparse.coo_matrix(np.asarray(vector, dtype=float).reshape(-1, 1))
    return canonical_sparse_fingerprints(
        matrix,
        names,
        ["value"],
        significant_digits,
    )


def _name_digest(names: list[str]) -> str:
    return hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()


def _deterministic_dual(name: str) -> float:
    token = int.from_bytes(hashlib.sha256(name.encode("utf-8")).digest()[:4], "big")
    return 0.5 + token / 2**32


def _design_variables(model: Any) -> dict[str, Any]:
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


def _parse_perturbation(specification: str) -> tuple[str, float]:
    try:
        name, delta = specification.split("=", 1)
        return name, float(delta)
    except ValueError as err:
        raise argparse.ArgumentTypeError(
            "Perturbations must have the form DESIGN_VARIABLE=DELTA."
        ) from err


def fingerprint(
    perturbations: list[tuple[str, float]],
    detail_output: Path | None = None,
) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    model.fs.c_tax_rate = 0.0
    checkpoint = _checkpoint(ARCHIVED_OPTIMUM)
    ms.from_json(model, fname=str(checkpoint))
    unfix_DOFs_pre_optimization(model)

    design = _design_variables(model)
    applied: list[dict[str, Any]] = []
    for name, delta in perturbations:
        if name not in design:
            raise ValueError(
                f"Unknown design variable {name!r}; choose from {sorted(design)}."
            )
        variable = design[name]
        before = float(value(variable))
        variable.set_value(before + delta)
        applied.append(
            {
                "key": name,
                "variable": variable.name,
                "delta": delta,
                "before": before,
                "after": float(value(variable)),
            }
        )

    jacobian, nlp = get_jacobian(model, scaled=False)
    variables = nlp.vlist
    constraints = nlp.clist
    variable_names = [variable.name for variable in variables]
    constraint_names = [constraint.name for constraint in constraints]
    gradient = nlp.evaluate_grad_objective()
    constraint_values = nlp.evaluate_constraints()

    duals = np.asarray([_deterministic_dual(name) for name in constraint_names])
    nlp.set_duals(duals)
    nlp.set_obj_factor(1.0)
    hessian = nlp.evaluate_hessian_lag()

    if detail_output is not None:
        detail_output.parent.mkdir(parents=True, exist_ok=True)
        jacobian_detail = sparse.csr_matrix(jacobian)
        jacobian_detail.eliminate_zeros()
        hessian_detail = sparse.csr_matrix(hessian)
        hessian_detail.eliminate_zeros()
        np.savez_compressed(
            detail_output,
            variable_names=np.asarray(variable_names),
            constraint_names=np.asarray(constraint_names),
            objective_gradient=np.asarray(gradient),
            constraint_values=np.asarray(constraint_values),
            jacobian_data=jacobian_detail.data,
            jacobian_indices=jacobian_detail.indices,
            jacobian_indptr=jacobian_detail.indptr,
            jacobian_shape=np.asarray(jacobian_detail.shape),
            hessian_data=hessian_detail.data,
            hessian_indices=hessian_detail.indices,
            hessian_indptr=hessian_detail.indptr,
            hessian_shape=np.asarray(hessian_detail.shape),
        )

    all_finite = all(
        math.isfinite(float(entry))
        for values in (jacobian.data, gradient, constraint_values, hessian.data)
        for entry in values
    )
    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-DERIVATIVE-FINGERPRINT-001",
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
            "perturbations": applied,
        },
        "all_derivative_and_constraint_values_finite": all_finite,
        "ordering": {
            "variable_names_sha256": _name_digest(variable_names),
            "constraint_names_sha256": _name_digest(constraint_names),
            "sorted_variable_names_sha256": _name_digest(sorted(variable_names)),
            "sorted_constraint_names_sha256": _name_digest(sorted(constraint_names)),
        },
        "detail_output": str(detail_output) if detail_output is not None else None,
        "objective_gradient": canonical_vector_fingerprints(
            gradient, variable_names
        ),
        "constraint_values": canonical_vector_fingerprints(
            constraint_values, constraint_names
        ),
        "constraint_jacobian": canonical_sparse_fingerprints(
            jacobian, constraint_names, variable_names
        ),
        "lagrangian_hessian": canonical_sparse_fingerprints(
            hessian, variable_names, variable_names
        ),
        "lagrangian_hessian_duals": {
            "definition": (
                "0.5 + uint32(big-endian, sha256(constraint_name)[:4]) / 2**32"
            ),
            "minimum": float(np.min(duals)),
            "maximum": float(np.max(duals)),
            "objective_factor": 1.0,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--perturb",
        action="append",
        default=[],
        type=_parse_perturbation,
        metavar="DESIGN_VARIABLE=DELTA",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--detail-output",
        type=Path,
        help="Optional compressed NumPy archive for name-aligned comparisons.",
    )
    args = parser.parse_args()
    report = fingerprint(args.perturb, detail_output=args.detail_output)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    _write_report(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
