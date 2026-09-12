#!/usr/bin/env python3
"""Compare actual row/column symbol ordering from two NL exports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def compare_sequences(
    reference: list[str], candidate: list[str], top: int = 20
) -> dict[str, Any]:
    reference_index = {name: index for index, name in enumerate(reference)}
    candidate_index = {name: index for index, name in enumerate(candidate)}
    shared = sorted(reference_index.keys() & candidate_index.keys())
    displacement = sorted(
        (
            abs(candidate_index[name] - reference_index[name]),
            name,
            reference_index[name],
            candidate_index[name],
        )
        for name in shared
    )
    prefix = 0
    for first, second in zip(reference, candidate):
        if first != second:
            break
        prefix += 1
    same_positions = sum(
        reference_index[name] == candidate_index[name] for name in shared
    )
    return {
        "reference_count": len(reference),
        "candidate_count": len(candidate),
        "component_sets_identical": set(reference) == set(candidate),
        "shared_components": len(shared),
        "same_position_count": same_positions,
        "longest_common_prefix": prefix,
        "maximum_absolute_position_displacement": (
            displacement[-1][0] if displacement else 0
        ),
        "largest_position_displacements": [
            {
                "name": name,
                "reference_index": first,
                "candidate_index": second,
                "absolute_displacement": delta,
            }
            for delta, name, first, second in reversed(displacement[-top:])
        ],
    }


def compare(reference_path: Path, candidate_path: Path, top: int = 20) -> dict[str, Any]:
    reference = json.loads(reference_path.read_text())
    candidate = json.loads(candidate_path.read_text())
    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-NL-SYMBOL-MAP-COMPARISON-001",
        "reference": {
            "path": str(reference_path),
            "environment": reference["environment"],
            "nl_file": reference["nl_file"],
        },
        "candidate": {
            "path": str(candidate_path),
            "environment": candidate["environment"],
            "nl_file": candidate["nl_file"],
        },
        "ordering": {
            category: compare_sequences(
                reference["ordering"][category]["names"],
                candidate["ordering"][category]["names"],
                top,
            )
            for category in ("variables", "constraints", "objectives")
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = compare(arguments.reference, arguments.candidate, arguments.top)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
