#!/usr/bin/env python3
"""Compare fresh Bakken ROK/tax figure payloads with workbook and paper."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from compare_regional_figure_data import (
    _compare_curves,
    _compare_frames,
    _curves_from_heat_table,
    _frame_from_payload,
)
from regenerate_archived_figures import C4_PLUS
from regenerate_composite_curves import DEFAULT_WORKBOOK, discover_cases


REPO_ROOT = Path(__file__).resolve().parents[1]


def _case_sheet(
    cases: list[Any], model_code: int, tax_rate: float
) -> Any:
    matching = [
        case for case in cases
        if case.model_code == model_code
        and case.region == "Bakken"
        and abs(case.tax_usd_per_kg - tax_rate) < 1e-12
    ]
    if len(matching) != 1:
        raise ValueError(
            f"Expected one Bakken M{model_code} tax={tax_rate} heat sheet; "
            f"found {len(matching)}."
        )
    return matching[0]


def _run_tax_rate(record: dict[str, Any], run: dict[str, Any]) -> float:
    tax_rate = run.get("co2_tax_usd_per_kg")
    if tax_rate is None:
        tax_rate = record["case"]["co2_tax_usd_per_kg"]
    return float(tax_rate)


def compare_records(
    records: list[dict[str, Any]],
    workbook_path: Path,
    reference: dict[str, Any],
) -> dict[str, Any]:
    cases = discover_cases(workbook_path)
    comparisons = []
    for record in records:
        model_code = record["case"]["model_code"]
        for run in record.get("runs") or [record.get("run", {})]:
            figure_data = run.get("figure_data")
            if figure_data is None:
                continue
            tax_rate = _run_tax_rate(record, run)
            case = _case_sheet(cases, model_code, tax_rate)
            stream_sheet = (
                f"streams_M{model_code}_tax={case.tax_token}_Bakken"
            )
            archived_streams = pd.read_excel(
                workbook_path, sheet_name=stream_sheet, index_col=0
            )
            archived_heat = pd.read_excel(
                workbook_path, sheet_name=case.sheet, index_col=0
            )
            fresh = run["results_in_migrated_csv_units"]
            archived_qw_mw = fresh["Qw"] - run[
                "difference_from_migrated_csv"
            ]["Qw"]
            archived_curves = _curves_from_heat_table(
                archived_heat, archived_qw_mw * 3.6
            )
            item: dict[str, Any] = {
                "model_code": model_code,
                "tax_usd_per_tonne": tax_rate * 1000.0,
                "termination_condition": run["termination_condition"],
                "stream_sheet": stream_sheet,
                "heat_integration_sheet": case.sheet,
                "stream_table": _compare_frames(
                    _frame_from_payload(figure_data["stream_table"]),
                    archived_streams,
                ),
                "heat_exchanger_table": _compare_frames(
                    _frame_from_payload(figure_data["heat_exchanger_table"]),
                    archived_heat,
                ),
                "composite_curves": _compare_curves(
                    figure_data["composite_curves"], archived_curves
                ),
            }
            if tax_rate == 0.045:
                model_name = f"M{model_code}"
                product = figure_data["liquid_product"]
                component_flow = product["component_flow_mol_per_s"]
                c4_plus_percent = (
                    100.0
                    * sum(component_flow.get(name, 0.0) for name in C4_PLUS)
                    / product["total_flow_mol_per_s"]
                )
                item["figure_s1"] = {
                    "fresh_total_mol_per_s": product["total_flow_mol_per_s"],
                    "paper_label_mol_per_s": reference["figure_s1"][
                        "total_labels"
                    ][model_name],
                    "difference_mol_per_s": (
                        product["total_flow_mol_per_s"]
                        - reference["figure_s1"]["total_labels"][model_name]
                    ),
                }
                item["figure_s2"] = {
                    "fresh_c4_plus_olefin_mole_percent": c4_plus_percent,
                    "paper_label_mole_percent": reference["figure_s2"][
                        "values"
                    ][model_name],
                    "difference_mole_percent": (
                        c4_plus_percent
                        - reference["figure_s2"]["values"][model_name]
                    ),
                }
                paper_emissions = next(
                    row for row in reference["figure_3"]["series"]
                    if row["case"] == model_name
                )
                item["figure_3"] = {
                    "upstream_difference_g_co2e_per_mj": (
                        figure_data["upstream_emissions_kg_co2e_per_gj"]
                        - paper_emissions["upstream"]
                    ),
                    "process_difference_g_co2e_per_mj": (
                        fresh["Downstream-em"] - paper_emissions["process"]
                    ),
                    "classification": "EMISSIONS-NORMALIZATION-001",
                }
            comparisons.append(item)
    return {
        "schema_version": 1,
        "workbook": str(workbook_path),
        "reference": "reproducibility/published_reference.json",
        "comparisons": comparisons,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-record", type=Path, action="append", required=True)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument(
        "--reference",
        type=Path,
        default=REPO_ROOT / "reproducibility" / "published_reference.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    report = compare_records(
        [json.loads(path.read_text()) for path in arguments.run_record],
        arguments.workbook,
        json.loads(arguments.reference.read_text()),
    )
    arguments.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(f"Compared {len(report['comparisons'])} fresh Bakken cases.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
