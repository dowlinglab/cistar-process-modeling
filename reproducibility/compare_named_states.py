#!/usr/bin/env python3
"""Compare two complete name-aligned variable-state snapshots."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any


def compare_values(reference: float, candidate: float) -> dict[str, float]:
    absolute = abs(candidate - reference)
    scale = max(1.0, abs(reference), abs(candidate))
    return {
        "absolute_difference": absolute,
        "scale_aware_difference": absolute / scale,
    }


def _component_family(name: str) -> str:
    return re.sub(r"\[[^]]*\]", "", name)


def _unit_name(name: str) -> str:
    fields = name.split(".")
    return ".".join(fields[:2]) if len(fields) >= 2 else name


def compare(
    reference_path: Path,
    candidate_path: Path,
    top: int,
    active_variable_map: Path | None = None,
) -> dict[str, Any]:
    reference_payload = json.loads(reference_path.read_text())
    candidate_payload = json.loads(candidate_path.read_text())
    reference = {item["name"]: item for item in reference_payload["variables"]}
    candidate = {item["name"]: item for item in candidate_payload["variables"]}
    selected_names = reference.keys() & candidate.keys()
    active_selection = None
    if active_variable_map is not None:
        active_payload = json.loads(active_variable_map.read_text())
        active_selection = set(active_payload["ordering"]["variables"]["names"])
        selected_names &= active_selection
    shared = sorted(selected_names)
    ranked = []
    nonfinite = []
    fixed_status_changes = []
    for name in shared:
        first = reference[name]
        second = candidate[name]
        if first["fixed"] != second["fixed"]:
            fixed_status_changes.append(name)
        first_value = first["value"]
        second_value = second["value"]
        if (
            first_value is None
            or second_value is None
            or not math.isfinite(first_value)
            or not math.isfinite(second_value)
        ):
            nonfinite.append(name)
            continue
        metrics = compare_values(first_value, second_value)
        ranked.append(
            (
                metrics["scale_aware_difference"],
                metrics["absolute_difference"],
                name,
                first_value,
                second_value,
            )
        )
    ranked.sort()
    changed = [item for item in ranked if item[1] != 0.0]
    thresholds = (1e-12, 1e-9, 1e-6, 1e-3, 1e-1)
    threshold_counts = {
        format(threshold, ".0e"): sum(item[0] > threshold for item in ranked)
        for threshold in thresholds
    }
    material = [item for item in ranked if item[0] > 1e-3]
    material_families = Counter(_component_family(item[2]) for item in material)
    material_units = Counter(_unit_name(item[2]) for item in material)

    def render(item: tuple[float, float, str, float, float]) -> dict[str, Any]:
        scaled, absolute, name, first, second = item
        return {
            "name": name,
            "reference_value": first,
            "candidate_value": second,
            "absolute_difference": absolute,
            "scale_aware_difference": scaled,
        }

    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-NAMED-STATE-COMPARISON-001",
        "reference": {
            "path": str(reference_path),
            "stage": reference_payload["stage"],
            "environment": reference_payload["environment"],
        },
        "candidate": {
            "path": str(candidate_path),
            "stage": candidate_payload["stage"],
            "environment": candidate_payload["environment"],
        },
        "counts": {
            "reference_variables": len(reference),
            "candidate_variables": len(candidate),
            "shared_variables": len(shared),
            "component_sets_identical": set(reference) == set(candidate),
            "exactly_equal_values": len(ranked) - len(changed),
            "different_values": len(changed),
            "nonfinite_or_missing_values": len(nonfinite),
            "fixed_status_changes": len(fixed_status_changes),
            "selected_variables": len(shared),
        },
        "selection": {
            "active_variable_map": (
                str(active_variable_map) if active_variable_map is not None else None
            ),
            "active_map_variable_count": (
                len(active_selection) if active_selection is not None else None
            ),
            "all_shared_variables_selected": active_variable_map is None,
        },
        "maxima": {
            "absolute_difference": render(max(ranked, key=lambda item: item[1])),
            "scale_aware_difference": render(ranked[-1]),
        },
        "scale_aware_threshold_counts": threshold_counts,
        "components_above_1e-3_by_family": [
            {"family": name, "count": count}
            for name, count in material_families.most_common()
        ],
        "components_above_1e-3_by_unit": [
            {"unit": name, "count": count}
            for name, count in material_units.most_common()
        ],
        "largest_scale_aware_differences": [
            render(item) for item in reversed(ranked[-top:])
        ],
        "largest_absolute_differences": [
            render(item)
            for item in reversed(sorted(ranked, key=lambda item: item[1])[-top:])
        ],
        "nonfinite_or_missing_names": nonfinite,
        "fixed_status_change_names": fixed_status_changes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument(
        "--active-variable-map",
        type=Path,
        help="Restrict comparison to names in an exported NL symbol-map JSON.",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = compare(
        arguments.reference,
        arguments.candidate,
        arguments.top,
        arguments.active_variable_map,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
