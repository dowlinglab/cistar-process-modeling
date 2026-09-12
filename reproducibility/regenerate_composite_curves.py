#!/usr/bin/env python3
"""Regenerate every archived optimal composite curve without overwriting it."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.utility_minimization_1d import gen_curves


DEFAULT_WORKBOOK = REPO_ROOT / "results" / "solution_data.xlsx"
SHEET_PATTERN = re.compile(
    r"^HI_M(?P<model_code>\d+)_tax=(?P<tax_token>.+?)_"
    r"(?P<region>.+)_optimal$"
)


@dataclass(frozen=True)
class CurveCase:
    sheet: str
    model_code: int
    tax_token: str
    tax_usd_per_kg: float
    region: str

    @property
    def filename(self) -> str:
        return (
            f"composite_curve_M{self.model_code}_C-tax_{self.tax_token}_"
            f"region_{self.region}_optimal.pdf"
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_case(sheet: str) -> CurveCase | None:
    match = SHEET_PATTERN.match(sheet)
    if match is None:
        return None
    region = match.group("region")
    if region == "EF-Basn":
        region = "EF-Basin"
    tax_token = match.group("tax_token")
    return CurveCase(
        sheet=sheet,
        model_code=int(match.group("model_code")),
        tax_token=tax_token,
        tax_usd_per_kg=float(tax_token),
        region=region,
    )


def discover_cases(workbook: Path) -> list[CurveCase]:
    return [
        case
        for sheet in pd.ExcelFile(workbook).sheet_names
        if (case := parse_case(sheet)) is not None
    ]


def load_result_tables(repo_root: Path = REPO_ROOT) -> dict[str, pd.DataFrame]:
    return {
        "models": pd.read_csv(repo_root / "optimal_data_wrt_ROK_models.csv"),
        "taxes": pd.read_csv(repo_root / "optimal_data_wrt_c_tax_rates.csv"),
        "regions": pd.read_csv(repo_root / "optimal_data_wrt_region.csv"),
    }


def cooling_utility_mw(
    case: CurveCase, tables: dict[str, pd.DataFrame]
) -> float:
    if case.region != "Bakken":
        table = tables["regions"]
        rows = table[(table["ROK_model"] == f"M{case.model_code}")]
        rows = rows[rows["Region"] == case.region]
    elif case.model_code != 5:
        table = tables["models"]
        rows = table[table["ROK_model"] == f"M{case.model_code}"]
    else:
        table = tables["taxes"]
        target = case.tax_usd_per_kg * 1000.0
        rows = table[np.isclose(table["C-tax-rate"], target, atol=1e-9)]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one cooling-utility row for {case.sheet}; found {len(rows)}"
        )
    return float(rows.iloc[0]["Qw"])


def curve_coordinates(
    frame: pd.DataFrame, cooling_utility: float
) -> dict[str, list[float]]:
    tin = pd.to_numeric(frame.iloc[0, 3:], errors="raise").to_numpy(float)
    tout = pd.to_numeric(frame.iloc[1, 3:], errors="raise").to_numpy(float)
    duty_w = pd.to_numeric(frame.iloc[2, 3:], errors="raise").to_numpy(float)
    duty_gj_per_hour = duty_w * 3.6e-6

    isothermal = np.isclose(tin, tout, rtol=0.0, atol=1e-9) & ~np.isclose(
        duty_gj_per_hour, 0.0
    )
    heating = (duty_gj_per_hour > 0) & ~isothermal
    cooling = (duty_gj_per_hour < 0) & ~isothermal
    t_hot, q_hot = gen_curves(
        tin[cooling], tout[cooling], duty_gj_per_hour[cooling]
    )
    t_cold, q_cold = gen_curves(
        tin[heating], tout[heating], -duty_gj_per_hour[heating]
    )
    q_cold = q_cold + duty_gj_per_hour[heating].sum()
    hot_mw = -q_hot * 1000.0 / 3600.0
    cold_mw = (q_cold + cooling_utility * 3.6) * 1000.0 / 3600.0

    hot_events = [
        (float(tin[i]), float(-duty_gj_per_hour[i] * 1000.0 / 3600.0))
        for i in np.flatnonzero(isothermal & (duty_gj_per_hour < 0))
    ]
    cold_events = [
        (float(tin[i]), float(duty_gj_per_hour[i] * 1000.0 / 3600.0))
        for i in np.flatnonzero(isothermal & (duty_gj_per_hour > 0))
    ]
    t_hot, hot_mw = add_isothermal_duties(t_hot, hot_mw, hot_events)
    t_cold, cold_mw = add_isothermal_duties(t_cold, cold_mw, cold_events)
    return {
        "hot_mw": hot_mw.tolist(),
        "hot_temperature_k": t_hot.tolist(),
        "cold_mw": cold_mw.tolist(),
        "cold_temperature_k": t_cold.tolist(),
    }


def add_isothermal_duties(
    temperatures: np.ndarray,
    cumulative_heat_mw: np.ndarray,
    events: list[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray]:
    """Insert zero-temperature-width heat duties as horizontal curve segments."""
    result_t = np.asarray(temperatures, dtype=float)
    result_q = np.asarray(cumulative_heat_mw, dtype=float)
    combined: dict[float, float] = {}
    for temperature, duty in events:
        combined[temperature] = combined.get(temperature, 0.0) + duty
    for temperature, duty in sorted(combined.items()):
        at_temperature = float(np.interp(temperature, result_t, result_q))
        below = result_t < temperature
        above = result_t > temperature
        result_t = np.concatenate(
            [result_t[below], [temperature, temperature], result_t[above]]
        )
        result_q = np.concatenate(
            [
                result_q[below],
                [at_temperature, at_temperature + duty],
                result_q[above] + duty,
            ]
        )
    return result_t, result_q


def save_plot(coordinates: dict[str, list[float]], output: Path) -> None:
    fig, axis = plt.subplots(figsize=(8, 6))
    axis.plot(
        coordinates["hot_mw"],
        coordinates["hot_temperature_k"],
        color="r",
        label="Hot Streams",
        linewidth=3,
    )
    axis.plot(
        coordinates["cold_mw"],
        coordinates["cold_temperature_k"],
        color="b",
        label="Cold Streams",
        linewidth=3,
    )
    axis.set_xlabel(
        "Cumulative process-wide \n heat exchange [MW]",
        fontsize=24,
        weight="bold",
    )
    axis.set_ylabel("Temperature [K]", fontsize=24, weight="bold")
    axis.tick_params(axis="both", labelsize=16)
    axis.legend(loc="best", fontsize=20)
    axis.grid()
    fig.savefig(output, bbox_inches="tight", dpi=200)
    plt.close(fig)


def regenerate(
    workbook: Path, output_dir: Path, record_path: Path | None
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = load_result_tables()
    report: dict[str, Any] = {
        "schema_version": 1,
        "source_workbook": str(workbook),
        "source_workbook_sha256": _sha256(workbook),
        "curve_count": 0,
        "curves": [],
    }
    for case in discover_cases(workbook):
        frame = pd.read_excel(workbook, sheet_name=case.sheet)
        qw_mw = cooling_utility_mw(case, tables)
        coordinates = curve_coordinates(frame, qw_mw)
        isothermal_duty_count = int(
            np.count_nonzero(
                np.isclose(
                    pd.to_numeric(frame.iloc[0, 3:]).to_numpy(float),
                    pd.to_numeric(frame.iloc[1, 3:]).to_numpy(float),
                    rtol=0.0,
                    atol=1e-9,
                )
                & ~np.isclose(pd.to_numeric(frame.iloc[2, 3:]).to_numpy(float), 0.0)
            )
        )
        output = output_dir / case.filename
        save_plot(coordinates, output)
        report["curves"].append(
            {
                "sheet": case.sheet,
                "model_code": case.model_code,
                "tax_usd_per_kg_co2e": case.tax_usd_per_kg,
                "region": case.region,
                "cooling_utility_mw": qw_mw,
                "isothermal_duty_count": isothermal_duty_count,
                "output": str(output),
                "coordinates": coordinates,
            }
        )
    report["curve_count"] = len(report["curves"])
    if record_path is not None:
        record_path.parent.mkdir(parents=True, exist_ok=True)
        record_path.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    report = regenerate(args.workbook.resolve(), args.output_dir.resolve(), args.record)
    print(
        json.dumps(
            {
                "curve_count": report["curve_count"],
                "output_dir": str(args.output_dir.resolve()),
                "source_workbook_sha256": report["source_workbook_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
