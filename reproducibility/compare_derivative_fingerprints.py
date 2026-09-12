#!/usr/bin/env python3
"""Compare two detailed derivative fingerprints after aligning component names."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse


def _load_matrix(data: Any, prefix: str) -> sparse.csr_matrix:
    return sparse.csr_matrix(
        (
            data[f"{prefix}_data"],
            data[f"{prefix}_indices"],
            data[f"{prefix}_indptr"],
        ),
        shape=tuple(data[f"{prefix}_shape"]),
    )


def _vector_entries(values: np.ndarray, names: np.ndarray) -> dict[tuple[str, ...], float]:
    return {
        (str(name),): float(entry)
        for name, entry in zip(names, values)
        if float(entry) != 0.0
    }


def _matrix_entries(
    matrix: sparse.spmatrix,
    row_names: np.ndarray,
    column_names: np.ndarray,
    symmetric: bool = False,
) -> dict[tuple[str, ...], float]:
    coordinate = matrix.tocoo()
    result: dict[tuple[str, ...], float] = {}
    for row, column, entry in zip(coordinate.row, coordinate.col, coordinate.data):
        row_name = str(row_names[row])
        column_name = str(column_names[column])
        key = (
            tuple(sorted((row_name, column_name)))
            if symmetric
            else (row_name, column_name)
        )
        result[key] = result.get(key, 0.0) + float(entry)
    return {key: entry for key, entry in result.items() if entry != 0.0}


def _compare(
    reference: dict[tuple[str, ...], float],
    candidate: dict[tuple[str, ...], float],
    top: int,
) -> dict[str, Any]:
    keys = sorted(reference.keys() | candidate.keys())
    differences = []
    for key in keys:
        first = reference.get(key, 0.0)
        second = candidate.get(key, 0.0)
        absolute = abs(second - first)
        relative = absolute / max(1.0, abs(first), abs(second))
        differences.append((relative, absolute, key, first, second))
    differences.sort(reverse=True)
    return {
        "reference_nonzeros": len(reference),
        "candidate_nonzeros": len(candidate),
        "union_entries": len(keys),
        "same_sparsity_pattern": set(reference) == set(candidate),
        "maximum_absolute_difference": max((item[1] for item in differences), default=0.0),
        "maximum_scaled_difference": max((item[0] for item in differences), default=0.0),
        "counts_above_scaled_tolerance": {
            str(tolerance): sum(item[0] > tolerance for item in differences)
            for tolerance in (1e-12, 1e-10, 1e-8, 1e-6)
        },
        "largest_scaled_differences": [
            {
                "components": list(key),
                "reference": first,
                "candidate": second,
                "absolute_difference": absolute,
                "scaled_difference": relative,
            }
            for relative, absolute, key, first, second in differences[:top]
        ],
    }


def compare(reference_path: Path, candidate_path: Path, top: int = 20) -> dict[str, Any]:
    with np.load(reference_path) as reference, np.load(candidate_path) as candidate:
        reference_variables = reference["variable_names"]
        candidate_variables = candidate["variable_names"]
        reference_constraints = reference["constraint_names"]
        candidate_constraints = candidate["constraint_names"]
        return {
            "schema_version": 1,
            "reference": str(reference_path),
            "candidate": str(candidate_path),
            "component_sets": {
                "variables_identical": set(reference_variables) == set(candidate_variables),
                "constraints_identical": set(reference_constraints) == set(candidate_constraints),
                "reference_variables": len(reference_variables),
                "candidate_variables": len(candidate_variables),
                "reference_constraints": len(reference_constraints),
                "candidate_constraints": len(candidate_constraints),
            },
            "objective_gradient": _compare(
                _vector_entries(reference["objective_gradient"], reference_variables),
                _vector_entries(candidate["objective_gradient"], candidate_variables),
                top,
            ),
            "constraint_values": _compare(
                _vector_entries(reference["constraint_values"], reference_constraints),
                _vector_entries(candidate["constraint_values"], candidate_constraints),
                top,
            ),
            "constraint_jacobian": _compare(
                _matrix_entries(
                    _load_matrix(reference, "jacobian"),
                    reference_constraints,
                    reference_variables,
                ),
                _matrix_entries(
                    _load_matrix(candidate, "jacobian"),
                    candidate_constraints,
                    candidate_variables,
                ),
                top,
            ),
            "lagrangian_hessian": _compare(
                _matrix_entries(
                    _load_matrix(reference, "hessian"),
                    reference_variables,
                    reference_variables,
                    symmetric=True,
                ),
                _matrix_entries(
                    _load_matrix(candidate, "hessian"),
                    candidate_variables,
                    candidate_variables,
                    symmetric=True,
                ),
                top,
            ),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = compare(arguments.reference, arguments.candidate, arguments.top)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
