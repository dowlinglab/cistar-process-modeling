#!/usr/bin/env python3
"""Compare fresh regional figure data with the archived workbook and paper."""

from __future__ import annotations

import argparse
import json
import numbers
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _frame_from_payload(payload: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        payload["data"],
        index=payload["index"],
        columns=payload["columns"],
    )


def _is_number(value: Any) -> bool:
    return isinstance(value, numbers.Real) and not pd.isna(value)


def _normalized_text(value: Any) -> str | None:
    if pd.isna(value):
        return None
    return str(value).strip()


def _compare_frames(fresh: pd.DataFrame, archived: pd.DataFrame) -> dict[str, Any]:
    fresh.index = fresh.index.map(str)
    fresh.columns = fresh.columns.map(str)
    archived.index = archived.index.map(str)
    archived.columns = archived.columns.map(str)

    shared_rows = sorted(set(fresh.index) & set(archived.index))
    shared_columns = sorted(set(fresh.columns) & set(archived.columns))
    numeric_differences = []
    text_mismatches = []
    unit_label_mismatches = []
    for row in shared_rows:
        for column in shared_columns:
            fresh_value = fresh.at[row, column]
            archived_value = archived.at[row, column]
            if _is_number(fresh_value) and _is_number(archived_value):
                absolute = abs(float(fresh_value) - float(archived_value))
                relative = absolute / max(abs(float(archived_value)), 1e-12)
                numeric_differences.append(
                    {
                        "row": row,
                        "column": column,
                        "fresh": float(fresh_value),
                        "archived": float(archived_value),
                        "absolute_difference": absolute,
                        "relative_difference": relative,
                    }
                )
            elif _normalized_text(fresh_value) != _normalized_text(archived_value):
                mismatch = {
                    "row": row,
                    "column": column,
                    "fresh": _normalized_text(fresh_value),
                    "archived": _normalized_text(archived_value),
                }
                if column in {"Unit", "Units"}:
                    unit_label_mismatches.append(mismatch)
                else:
                    text_mismatches.append(mismatch)

    by_absolute = sorted(
        numeric_differences,
        key=lambda item: item["absolute_difference"],
        reverse=True,
    )
    by_relative = sorted(
        numeric_differences,
        key=lambda item: item["relative_difference"],
        reverse=True,
    )
    material_relative_floor = 1e-8
    material_by_relative = [
        item
        for item in by_relative
        if abs(item["archived"]) >= material_relative_floor
    ]
    return {
        "fresh_shape": list(fresh.shape),
        "archived_shape": list(archived.shape),
        "missing_rows_in_fresh": sorted(set(archived.index) - set(fresh.index)),
        "extra_rows_in_fresh": sorted(set(fresh.index) - set(archived.index)),
        "missing_columns_in_fresh": sorted(
            set(archived.columns) - set(fresh.columns)
        ),
        "extra_columns_in_fresh": sorted(
            set(fresh.columns) - set(archived.columns)
        ),
        "numeric_cells_compared": len(numeric_differences),
        "maximum_absolute_difference": (
            by_absolute[0]["absolute_difference"] if by_absolute else None
        ),
        "maximum_relative_difference": (
            by_relative[0]["relative_difference"] if by_relative else None
        ),
        "material_relative_difference_floor": material_relative_floor,
        "maximum_material_relative_difference": (
            material_by_relative[0]["relative_difference"]
            if material_by_relative
            else None
        ),
        "largest_absolute_differences": by_absolute[:10],
        "largest_relative_differences": by_relative[:10],
        "text_mismatches": text_mismatches[:20],
        "text_mismatch_count": len(text_mismatches),
        "unit_label_mismatches": unit_label_mismatches[:20],
        "unit_label_mismatch_count": len(unit_label_mismatches),
    }


def _paper_lookup(reference: dict[str, Any], region: str) -> dict[str, float]:
    for row in reference["figure_6"]["series"]:
        if row["region"] == region:
            return {"upstream": row["upstream"], "process": row["process"]}
    raise KeyError(f"No Figure 6 row for {region}")


def _workbook_region(region: str) -> str:
    return "EF-Basn" if region == "EF-Basin" else region


def _gen_curves(tin: np.ndarray, tout: np.ndarray, heat: np.ndarray):
    temperatures = np.unique(np.concatenate((tin, tout)))
    cumulative = np.zeros(len(temperatures))
    for inlet, outlet, duty in zip(tin, tout, heat):
        x = [inlet, outlet]
        y = [0.0, duty]
        lower = 0 if inlet < outlet else 1
        upper = 1 - lower
        for index, temperature in enumerate(temperatures):
            alpha = (x[upper] - temperature) / (x[upper] - x[lower])
            alpha = min(max(alpha, 0.0), 1.0)
            cumulative[index] += alpha * (y[upper] - y[lower]) + y[lower]
    return temperatures, cumulative


def _curves_from_heat_table(
    heat_table: pd.DataFrame, cooling_utility_gj_per_hour: float
) -> dict[str, dict[str, list[float]]]:
    rows = heat_table.set_index("Quantity")
    heating = ["fs.H101", "fs.H103", "fs.R101"]
    cooling = ["fs.H102", "fs.H104", "fs.H105", "fs.H106", "fs.R102"]

    heating_inlet = rows.loc["T inlet", heating].astype(float).to_numpy()
    heating_outlet = (
        rows.loc["T outlet", heating].astype(float).to_numpy() + 1.0
    )
    heating_duty = (
        rows.loc["Heat duty", heating].astype(float).to_numpy() * 3.6e-6
    )
    cooling_inlet = (
        rows.loc["T inlet", cooling].astype(float).to_numpy() + 1.0
    )
    cooling_outlet = rows.loc["T outlet", cooling].astype(float).to_numpy()
    cooling_duty = (
        rows.loc["Heat duty", cooling].astype(float).to_numpy() * 3.6e-6
    )

    hot_temperature, hot_heat = _gen_curves(
        cooling_inlet, cooling_outlet, cooling_duty
    )
    cold_temperature, cold_heat = _gen_curves(
        heating_inlet, heating_outlet, -heating_duty
    )
    cold_heat = cold_heat + sum(heating_duty) + cooling_utility_gj_per_hour
    return {
        "hot": {
            "temperature_k": hot_temperature.tolist(),
            "cumulative_heat_gj_per_hour": (-hot_heat).tolist(),
        },
        "cold": {
            "temperature_k": cold_temperature.tolist(),
            "cumulative_heat_gj_per_hour": cold_heat.tolist(),
        },
    }


def _compare_curves(
    fresh: dict[str, dict[str, list[float]]],
    archived: dict[str, dict[str, list[float]]],
) -> dict[str, Any]:
    result = {}
    for side in ("hot", "cold"):
        fresh_temperature = np.asarray(fresh[side]["temperature_k"], dtype=float)
        archived_temperature = np.asarray(
            archived[side]["temperature_k"], dtype=float
        )
        fresh_heat = np.asarray(
            fresh[side]["cumulative_heat_gj_per_hour"], dtype=float
        )
        archived_heat = np.asarray(
            archived[side]["cumulative_heat_gj_per_hour"], dtype=float
        )
        comparison_grid = np.unique(
            np.concatenate((fresh_temperature, archived_temperature))
        )
        fresh_interpolated = np.interp(
            comparison_grid, fresh_temperature, fresh_heat
        )
        archived_interpolated = np.interp(
            comparison_grid, archived_temperature, archived_heat
        )
        nearest_knot_distances = [
            min(abs(point - archived_temperature))
            for point in fresh_temperature
        ] + [
            min(abs(point - fresh_temperature))
            for point in archived_temperature
        ]
        result[side] = {
            "fresh_point_count": len(fresh_temperature),
            "archived_point_count": len(archived_temperature),
            "comparison_grid_point_count": len(comparison_grid),
            "maximum_nearest_knot_distance_k": max(
                nearest_knot_distances, default=None
            ),
            "maximum_heat_difference_gj_per_hour": float(
                np.max(abs(fresh_interpolated - archived_interpolated))
            ),
        }
    return result


def _figure_7_classification(
    difference: float, stream_max_material_relative: float | None
) -> str:
    if abs(difference) <= 0.5:
        return "matches_printed_label"
    if (
        stream_max_material_relative is not None
        and stream_max_material_relative <= 1e-5
    ):
        return "PUBLISHED-SNAPSHOT-DRIFT-001"
    return "unclassified_fresh_solution_difference"


def compare_record(
    record: dict[str, Any], workbook_path: Path, reference: dict[str, Any]
) -> dict[str, Any]:
    source_runs = record.get("runs") or [record.get("run", {})]
    comparisons = []
    for source_run in source_runs:
        figure_data = source_run.get("figure_data")
        if figure_data is None:
            continue
        region = source_run.get("region") or record["case"]["region"]
        workbook_region = _workbook_region(region)
        stream_sheet = f"streams_M5_tax=0.045_{region}"
        heat_sheet = f"HI_M5_tax=0.045_{workbook_region}_optimal"
        archived_streams = pd.read_excel(
            workbook_path, sheet_name=stream_sheet, index_col=0
        )
        archived_heat = pd.read_excel(
            workbook_path, sheet_name=heat_sheet, index_col=0
        )

        paper_emissions = _paper_lookup(reference, region)
        fresh_upstream = figure_data["upstream_emissions_kg_co2e_per_gj"]
        fresh_process = source_run["results_in_migrated_csv_units"][
            "Downstream-em"
        ]
        fresh_lhv = sum(
            figure_data["liquid_product"][
                "component_lhv_contribution_mj_per_s"
            ].values()
        )
        paper_lhv = reference["figure_7"]["total_labels"][region]
        archived_qw_mw = (
            source_run["results_in_migrated_csv_units"]["Qw"]
            - source_run["difference_from_migrated_csv"]["Qw"]
        )
        archived_curves = _curves_from_heat_table(
            archived_heat, archived_qw_mw * 3600.0 / 1000.0
        )
        stream_comparison = _compare_frames(
            _frame_from_payload(figure_data["stream_table"]),
            archived_streams,
        )
        heat_comparison = _compare_frames(
            _frame_from_payload(figure_data["heat_exchanger_table"]),
            archived_heat,
        )
        lhv_difference = fresh_lhv - paper_lhv
        emissions_match_printed = (
            abs(fresh_upstream - paper_emissions["upstream"]) <= 0.005
            and abs(fresh_process - paper_emissions["process"]) <= 0.005
        )

        comparisons.append(
            {
                "region": region,
                "stream_sheet": stream_sheet,
                "heat_integration_sheet": heat_sheet,
                "stream_table": stream_comparison,
                "heat_exchanger_table": heat_comparison,
                "composite_curves": _compare_curves(
                    figure_data["composite_curves"], archived_curves
                ),
                "figure_6": {
                    "unit": "g CO2e/MJ fuel",
                    "fresh_upstream": fresh_upstream,
                    "paper_upstream": paper_emissions["upstream"],
                    "upstream_difference": (
                        fresh_upstream - paper_emissions["upstream"]
                    ),
                    "fresh_process": fresh_process,
                    "paper_process": paper_emissions["process"],
                    "process_difference": (
                        fresh_process - paper_emissions["process"]
                    ),
                    "within_print_rounding": emissions_match_printed,
                    "classification": (
                        "matches_printed_labels"
                        if emissions_match_printed
                        else "EMISSIONS-NORMALIZATION-001"
                    ),
                },
                "figure_7": {
                    "unit": "MW liquid hydrocarbon LHV",
                    "fresh_total": fresh_lhv,
                    "paper_total_label": paper_lhv,
                    "difference": lhv_difference,
                    "within_print_rounding": abs(lhv_difference) <= 0.5,
                    "classification": _figure_7_classification(
                        lhv_difference,
                        stream_comparison[
                            "maximum_material_relative_difference"
                        ],
                    ),
                },
            }
        )

    return {
        "schema_version": 1,
        "source_run_case": record.get("case"),
        "comparison_targets": {
            "workbook": str(workbook_path),
            "paper_reference": "reproducibility/published_reference.json",
        },
        "regions_with_figure_data": [item["region"] for item in comparisons],
        "comparisons": comparisons,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-record", type=Path, required=True)
    parser.add_argument(
        "--workbook",
        type=Path,
        default=REPO_ROOT / "results" / "solution_data.xlsx",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=REPO_ROOT / "reproducibility" / "published_reference.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    record = json.loads(args.run_record.read_text(encoding="utf-8"))
    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    report = compare_record(record, args.workbook, reference)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
