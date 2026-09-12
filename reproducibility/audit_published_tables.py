#!/usr/bin/env python3
"""Audit checked-in optimization summaries against published SI tables.

This audit intentionally keeps the printed paper values separate from the
higher-precision repository snapshots.  It reports known discrepancies rather
than changing either source to make the comparison pass.
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

CSV_FIELDS = {
    "msp_usd_per_mj": ("MSP", 1.0),
    "downstream_emissions_g_co2e_per_mj": ("Downstream-em", 1.0),
    "product_lhv_mw": ("Product-LHV", 1.0),
    "h2_rebate_usd_per_year": ("H2-rebate", 1.0),
    "tac_usd_per_year": ("TAC", 1.0),
    "t_r102_k": ("T_R102", 1.0),
    "t_h104_k": ("T_H104", 1.0),
    "t_h105_k": ("T_H105", 1.0),
    "p_f101_bar": ("P_F101", 1.0e-5),
    "t_h106_k": ("T_H106", 1.0),
    "p_h106_bar": ("P_H106", 1.0e-5),
    "p_f102_bar": ("P_F102", 1.0e-5),
    "qs_mw": ("Qs", 1.0),
    "qw_mw": ("Qw", 1.0),
}

SNAPSHOT_FILES = {
    "postprocessed": {
        "table_s4": REPO_ROOT / "optimal_data_wrt_ROK_models.csv",
        "table_s5": REPO_ROOT / "optimal_data_wrt_c_tax_rates.csv",
        "table_s6": REPO_ROOT / "optimal_data_wrt_region.csv",
    },
    "migrated": {
        "table_s4": REPO_ROOT / "results" / "optimal_data_wrt_ROK_models.csv",
        "table_s5": REPO_ROOT / "results" / "optimal_data_wrt_c_tax_rates.csv",
        "table_s6": REPO_ROOT / "results" / "optimal_data_wrt_region.csv",
    },
}

EMISSIONS_ISSUE = "EMISSIONS-NORMALIZATION-001"
PUBLISHED_SNAPSHOT_ISSUE = "PUBLISHED-SNAPSHOT-DRIFT-001"
LEGACY_ECONOMICS_ISSUE = "LEGACY-ECONOMICS-POSTPROCESSING-001"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _row_key(table_name: str, row: dict[str, Any]) -> str:
    if table_name == "table_s4":
        return str(row["rok_model"])
    if table_name == "table_s5":
        return format(float(row["co2_tax_usd_per_tonne"]), "g")
    if table_name == "table_s6":
        return str(row["region"])
    raise ValueError(f"Unsupported table: {table_name}")


def _csv_row_key(table_name: str, row: dict[str, str]) -> str:
    if table_name == "table_s4":
        return row["ROK_model"]
    if table_name == "table_s5":
        return format(float(row["C-tax-rate"]), "g")
    if table_name == "table_s6":
        return row["Region"]
    raise ValueError(f"Unsupported table: {table_name}")


def audit_snapshot(snapshot: str = "postprocessed") -> dict[str, Any]:
    """Return a structured comparison for one checked-in snapshot."""

    reference = _load_json(REFERENCE_PATH)
    tolerances = reference["comparison_policy"]["metrics"]
    if snapshot not in SNAPSHOT_FILES:
        raise ValueError(f"Unknown snapshot {snapshot!r}")

    comparisons: list[dict[str, Any]] = []
    missing_rows: list[dict[str, str]] = []

    for table_name, csv_path in SNAPSHOT_FILES[snapshot].items():
        published_rows = reference[table_name]["rows"]
        archived_rows = {
            _csv_row_key(table_name, row): row for row in _load_csv(csv_path)
        }

        for published in published_rows:
            key = _row_key(table_name, published)
            archived = archived_rows.get(key)
            if archived is None:
                missing_rows.append(
                    {"table": table_name, "row": key, "csv": str(csv_path)}
                )
                continue

            for metric, (csv_field, multiplier) in CSV_FIELDS.items():
                expected = float(published[metric])
                actual = float(archived[csv_field]) * multiplier
                difference = actual - expected
                tolerance = float(tolerances[metric])
                matches = math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance)
                known_issue = None
                if not matches:
                    if metric == "downstream_emissions_g_co2e_per_mj":
                        known_issue = EMISSIONS_ISSUE
                    elif snapshot == "migrated" and metric in {
                        "msp_usd_per_mj",
                        "tac_usd_per_year",
                    }:
                        known_issue = LEGACY_ECONOMICS_ISSUE
                    else:
                        known_issue = PUBLISHED_SNAPSHOT_ISSUE
                comparisons.append(
                    {
                        "table": table_name,
                        "row": key,
                        "metric": metric,
                        "published": expected,
                        "snapshot": actual,
                        "difference": difference,
                        "absolute_difference": abs(difference),
                        "tolerance": tolerance,
                        "matches": matches,
                        "known_issue": known_issue,
                    }
                )

    mismatches = [item for item in comparisons if not item["matches"]]
    unexpected = [item for item in mismatches if item["known_issue"] is None]
    known = [item for item in mismatches if item["known_issue"] is not None]
    return {
        "schema_version": 1,
        "snapshot": snapshot,
        "reference": str(REFERENCE_PATH.relative_to(REPO_ROOT)),
        "comparison_count": len(comparisons),
        "match_count": sum(item["matches"] for item in comparisons),
        "mismatch_count": len(mismatches),
        "known_mismatch_count": len(known),
        "unexpected_mismatch_count": len(unexpected),
        "missing_row_count": len(missing_rows),
        "known_issues": {
            EMISSIONS_ISSUE: {
                "status": "unresolved",
                "summary": (
                    "Published process/downstream emissions do not match the "
                    "checked-in CSV values. The public source normalizes by "
                    "fuel_energy + gas_energy + h2_energy; the paper values "
                    "appear to use an earlier convention absent from public "
                    "Git history. Do not change the reference until the exact "
                    "published denominator is recovered or reconstructed."
                ),
            },
            PUBLISHED_SNAPSHOT_ISSUE: {
                "status": "unresolved",
                "summary": (
                    "A small subset of printed SI values does not equal the "
                    "checked-in higher-precision result rounded to the displayed "
                    "precision. This is evidence that the publication tables and "
                    "public CSV snapshot are not the exact same optimization "
                    "snapshot; the public Git history begins with a bulk migration "
                    "and cannot identify the earlier run."
                ),
            },
            LEGACY_ECONOMICS_ISSUE: {
                "status": "explained",
                "summary": (
                    "The migrated results/ CSVs predate commit 957e363, which "
                    "recalculated TAC and MSP using separate cooling-water and "
                    "refrigerated-water grades. The root-level postprocessed CSVs "
                    "contain that publication-aligned economics correction."
                ),
            },
        },
        "missing_rows": missing_rows,
        "mismatches": mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapshot",
        choices=sorted(SNAPSHOT_FILES),
        default="postprocessed",
        help="Checked-in result family to audit.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a nonzero status for any mismatch, including known issues.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for the JSON report; stdout is always populated.",
    )
    args = parser.parse_args()

    report = audit_snapshot(args.snapshot)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")

    if report["missing_row_count"] or report["unexpected_mismatch_count"]:
        return 2
    if args.strict and report["mismatch_count"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
