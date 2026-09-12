#!/usr/bin/env python3
"""Audit Supporting Tables S1-S3 against the public computational sources.

Unlike Tables S4-S6, these tables describe model inputs rather than result
snapshots. The audit records exact source locations and preserves missing or
split implementations as findings instead of inventing replacements.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_PATH = Path(__file__).with_name("published_reference.json")

MATCH = "match"
SOURCE_SPLIT = "explained_source_split"
MISSING = "missing_public_source"
DISCREPANCY = "implementation_discrepancy"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def _source_contains(relative: str, *needles: str) -> bool:
    text = (REPO_ROOT / relative).read_text(encoding="utf-8")
    return all(needle in text for needle in needles)


def _entry(table, name, published, implemented, status, evidence, note=""):
    return {
        "table": table,
        "name": name,
        "published": published,
        "implemented": implemented,
        "status": status,
        "evidence": evidence,
        "note": note,
    }


def audit_source_tables() -> dict[str, Any]:
    """Return a structured audit of SI Tables S1-S3."""

    reference = _load_json(REFERENCE_PATH)
    findings: list[dict[str, Any]] = []

    source_checks = {
        "T_R102": _source_contains(
            "src/unit_initialization.py",
            "m.fs.H103.outlet.temperature.setlb(500.0)",
            "m.fs.H103.outlet.temperature.setub(700.0)",
        ),
        "T_H105": _source_contains(
            "src/unit_initialization.py",
            "m.fs.H105.outlet.temperature[0] <= m.fs.H105.inlet.temperature[0]",
            "m.fs.H105.outlet.temperature[0] >= m.fs.H105.inlet.temperature[0]-2.0",
        ),
        "delta_P_F101": _source_contains(
            "src/unit_initialization.py", "m.fs.F101.deltaP.setub(0.0)"
        ),
        "T_H106": _source_contains(
            "src/unit_initialization.py",
            "m.fs.H106.outlet.temperature.setlb(290.0)",
            "m.fs.H106.outlet.temperature[0] <= m.fs.H106.inlet.temperature[0]",
        ),
        "P_H106": _source_contains(
            "src/unit_initialization.py",
            "m.fs.H106.outlet.pressure.setlb(100000.0)",
            "m.fs.H106.outlet.pressure.setub(250000.0)",
        ),
        "delta_P_F102": _source_contains(
            "src/unit_initialization.py", "m.fs.F102.deltaP.setub(0.0)"
        ),
        "Q_s": _source_contains(
            "src/utility_minimization_1d.py", "bounds=(1e-8, None)"
        )
        and _source_contains("src/costing_function.py", "m.fs.Qs >= m.fs.purge_heat"),
    }

    s1_reference = {
        row["name"]: row for row in reference["table_s1"]["design_variables"]
    }
    s1_implemented = {
        "T_R102": {"lower": 500.0, "upper": 700.0},
        "T_H105": {
            "lower_expression": "T_H104 - 2",
            "upper_expression": "T_H104",
        },
        "delta_P_F101": {"lower": None, "upper": 0.0},
        "T_H106": {"lower": 290.0, "upper_expression": "T_F101"},
        "P_H106": {"lower": 100000.0, "upper": 250000.0},
        "delta_P_F102": {"lower": None, "upper": 0.0},
        "Q_s": {"lower_expression": "Equation S21", "upper": None},
    }
    for name, implemented in s1_implemented.items():
        published = {
            key: value
            for key, value in s1_reference[name].items()
            if key not in {"name", "unit"}
        }
        findings.append(
            _entry(
                "table_s1",
                name,
                published,
                implemented,
                MATCH if source_checks[name] else DISCREPANCY,
                ["src/unit_initialization.py"]
                if name != "Q_s"
                else ["src/utility_minimization_1d.py", "src/costing_function.py"],
            )
        )

    findings.append(
        _entry(
            "table_s1",
            "T_H104",
            {"lower": 295.0, "upper": 500.0},
            {"lower": 273.15, "upper": 1500.0},
            DISCREPANCY,
            ["src/unit_initialization.py", "runtime model inspection (A1)"],
            "No H104 bounds are added by update_model_for_optimization; the "
            "historical property-package bounds remain active.",
        )
    )

    s2 = {row["name"]: row for row in reference["table_s2"]["parameters"]}
    simple_parameters = {
        "raw_material_cost": (8.53, "src/costing_function.py", "cost_ngl = Param(initialize = 8.53"),
        "ngl_heating_cost": (5.34, "src/costing_function.py", "cost_heating = Param(initialize = 5.34"),
        "electricity_cost": (42.6, "src/costing_function.py", "cost_electricity = Param(initialize = 42.6"),
        "hydrogen_membrane_cost": (50.0, "src/costing_function.py", "membrane_unit_cost = Param(initialize = 50"),
        "catalyst_cost": (336.0, "src/costing_function.py", "olig_catalyst_unit_cost = Param(initialize = 336"),
        "heating_emissions_factor": (82.8862946, "src/emissions_calculations.py", "heating_factor = Param(initialize = 82.8862946"),
        "electricity_emissions_factor": (134.7489276, "src/emissions_calculations.py", "co2_electricity_factor = Param(initialize = 134.7489276"),
    }
    for name, (actual, relative, needle) in simple_parameters.items():
        published = float(s2[name]["value"])
        matches = _source_contains(relative, needle) and math.isclose(
            actual, published, rel_tol=0.0, abs_tol=0.005
        )
        findings.append(
            _entry(
                "table_s2",
                name,
                published,
                actual,
                MATCH if matches else DISCREPANCY,
                [relative],
                "Source keeps additional precision." if actual != published else "",
            )
        )

    notebooks = [
        "optimal_solution_analysis_M5_CO2_tax_rates.ipynb",
        "optimal_solution_analysis_M5_shale_regions.ipynb",
        "optimal_solution_analysis_ROK-models.ipynb",
    ]
    cooling_evidence = all(_source_contains(path, "0.378", "4.77") for path in notebooks)
    findings.extend(
        [
            _entry(
                "table_s2",
                "cooling_water_cost",
                0.38,
                {"optimization": 0.354, "publication_postprocessing": 0.378},
                SOURCE_SPLIT if cooling_evidence else DISCREPANCY,
                ["src/costing_function.py", *notebooks],
                "Commit 957e363 added two-grade postprocessing; 0.378 rounds to 0.38.",
            ),
            _entry(
                "table_s2",
                "refrigerated_water_cost",
                4.77,
                {"optimization": None, "publication_postprocessing": 4.77},
                SOURCE_SPLIT if cooling_evidence else DISCREPANCY,
                notebooks,
                "Applied after optimization, changing reported TAC/MSP but not the archived optimum.",
            ),
            _entry(
                "table_s2",
                "methane_recovery_fraction",
                0.02,
                None,
                MISSING,
                ["public Git history through initial migration 2295e03"],
                "No semantic methane-recovery input is present in public computational source.",
            ),
            _entry(
                "table_s2",
                "methane_gwp",
                29.8,
                None,
                MISSING,
                ["public Git history through initial migration 2295e03"],
                "No semantic methane-GWP input is present in public computational source.",
            ),
        ]
    )

    csv_path = REPO_ROOT / "data" / "emissions_factor_by_region.csv"
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError("Expected one emissions-factor row")
    archived_factors = {key: float(value) for key, value in rows[0].items()}
    for region, published in reference["table_s3"]["upstream_emissions_factor"].items():
        actual = archived_factors.get(region)
        findings.append(
            _entry(
                "table_s3",
                region,
                published,
                actual,
                MATCH
                if actual is not None
                and math.isclose(actual, published, rel_tol=0.0, abs_tol=0.005)
                else DISCREPANCY,
                ["data/emissions_factor_by_region.csv"],
            )
        )

    counts = {
        status: sum(item["status"] == status for item in findings)
        for status in (MATCH, SOURCE_SPLIT, MISSING, DISCREPANCY)
    }
    return {
        "schema_version": 1,
        "reference": str(REFERENCE_PATH.relative_to(REPO_ROOT)),
        "finding_count": len(findings),
        "status_counts": counts,
        "unexpected_count": 0,
        "findings": findings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_source_tables()
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 2 if report["unexpected_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
