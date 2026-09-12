#!/usr/bin/env python3
"""Regenerate the eight remaining quantitative figures from archived sources.

The saved optimal checkpoints are authoritative for model-derived composition,
LHV, and emissions series. Checked-in CSV inputs are authoritative for the two
MSP plots and feed-composition plot. Outputs must be written outside the repo;
the publication PDFs are never overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from idaes.core.util import model_serializer as ms
from pyomo.environ import value

from reproducibility.run_m5_bakken_tax_series import _build_preoptimization_model
from src.plotting_functions import colors_dict


REFERENCE_PATH = Path(__file__).with_name("published_reference.json")
ROK_MODELS = ("M2", "M3", "M4", "M5")
REGIONS = (
    "Bakken",
    "EF-1",
    "EF-2",
    "EF-3",
    "EF-4",
    "EF-5",
    "EF-6",
    "EF-7",
    "EF-8",
    "EF-9",
    "EF-10",
    "EF-11",
    "EF-12",
    "EF-Basin",
)
C4_PLUS = ("butene", "pentene", "hexene", "heptene", "octene", "nonene")
OUTPUT_NAMES = {
    "F3": "emissions_by_ROK_models.pdf",
    "F5": "msp_vs_C-tax_rate_Bakken_M5.pdf",
    "F6": "emissions_by_region.pdf",
    "F7": "LHV_contribution_by_region.pdf",
    "F8": "msp_by_region_scatter_with_labels.pdf",
    "FS1": "outlet_component_flowrate_by_model.pdf",
    "FS2": "outlet_comp_C4-C9_alkenes_ROK_models.pdf",
    "FS5": "feed_compositions_bakken_eagleford.pdf",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rok_checkpoint(model_name: str) -> Path:
    suffix = "_sequential_solve" if model_name == "M5" else ""
    return REPO_ROOT / "initialization_files" / (
        f"CISTAR_optimal_solution_Bakken_C_tax_0.045_{model_name}_purge_0.01"
        f"{suffix}.json.gz"
    )


def region_checkpoint(region: str) -> Path:
    return REPO_ROOT / "initialization_files" / (
        f"CISTAR_optimal_solution_{region}_C_tax_0.045_M5_purge_0.01_"
        "sequential_solve.json.gz"
    )


def c4_plus_percent(flows: dict[str, float]) -> float:
    total = sum(flows.values())
    if total <= 0:
        raise ValueError("Liquid flow total must be positive")
    return 100.0 * sum(flows.get(component, 0.0) for component in C4_PLUS) / total


def _extract_model_state(model: Any) -> dict[str, Any]:
    components = list(model.fs.liquid)
    flows = OrderedDict(
        (
            component,
            value(model.fs.F102.liq_outlet.flow_mol_phase_comp[0, "Liq", component]),
        )
        for component in components
    )
    total = sum(flows.values())
    lhv = OrderedDict(
        (
            component,
            value(model.fs.LHV_per_component_per_stream["liq_outlet", component])
            * total,
        )
        for component in components
    )
    return {
        "liquid_component_flow_mol_per_s": flows,
        "total_liquid_flow_mol_per_s": total,
        "c4_plus_olefin_mole_percent": c4_plus_percent(flows),
        "liquid_component_lhv_mw": lhv,
        "total_liquid_lhv_mw": sum(lhv.values()),
        "upstream_emissions_g_co2e_per_mj": value(model.fs.upstream_emissions),
        "process_emissions_g_co2e_per_mj": value(model.fs.downstream_emissions),
    }


def collect_checkpoint_data() -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    """Build each topology once and load all archived optimum snapshots."""

    rok: dict[str, Any] = OrderedDict()
    hashes: dict[str, str] = OrderedDict()
    m5_model = None
    for code in (2, 3, 4, 5):
        model_name = f"M{code}"
        model = _build_preoptimization_model(
            model_code=code,
            region="Bakken",
            costing_tax=0.0 if code == 5 else 0.045,
        )
        checkpoint = rok_checkpoint(model_name)
        ms.from_json(model, fname=str(checkpoint))
        rok[model_name] = _extract_model_state(model)
        hashes[str(checkpoint.relative_to(REPO_ROOT))] = _sha256(checkpoint)
        if code == 5:
            m5_model = model

    assert m5_model is not None
    regional: dict[str, Any] = OrderedDict()
    for region in REGIONS:
        checkpoint = region_checkpoint(region)
        ms.from_json(m5_model, fname=str(checkpoint))
        regional[region] = _extract_model_state(m5_model)
        hashes[str(checkpoint.relative_to(REPO_ROOT))] = _sha256(checkpoint)
    return rok, regional, hashes


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _save(fig: Any, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", dpi=200)
    plt.close(fig)


def _stacked_emissions(
    labels: list[str], upstream: list[float], process: list[float], path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(labels))
    p1 = ax.bar(x, upstream, width=0.75, label="Upstream")
    p2 = ax.bar(x, process, width=0.75, bottom=upstream, label="Process")
    ax.set_ylabel("GHG emissions\n[g CO$_2$e/MJ fuel]", fontsize=24, weight="bold")
    ax.set_xlabel(
        "Oligomerization model" if len(labels) == 5 else "Shale region",
        fontsize=24,
        weight="bold",
    )
    ax.set_xticks(x, labels, fontsize=16, rotation=45)
    ax.tick_params(axis="y", labelsize=16)
    ax.bar_label(p1, label_type="center", fontsize=14, weight="bold", rotation=90, fmt="%.2f")
    ax.bar_label(p2, label_type="center", fontsize=14, weight="bold", rotation=90, fmt="%.2f")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.34), ncol=2, fontsize=20)
    _save(fig, path)


def _component_bars(
    data: dict[str, Any], labels: tuple[str, ...], key: str, horizontal: bool, path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    components = list(data[labels[0]][key])
    cumulative = np.zeros(len(labels))
    for component in components:
        values = [data[label][key][component] for label in labels]
        if horizontal:
            bars = ax.barh(labels, values, 0.5, left=cumulative, label=component, color=colors_dict[component])
        else:
            bars = ax.bar(labels, values, 0.5, bottom=cumulative, label=component, color=colors_dict[component])
        cumulative += values
    if horizontal:
        ax.set_xlabel("Liquid hydrocarbon flow rate [mol/s]", fontsize=24, weight="bold")
        ax.set_ylabel("ROK model", fontsize=24, weight="bold")
        ax.set_xlim(0, 210)
        ax.bar_label(bars, label_type="edge", fontsize=16, weight="bold", fmt="%.1f")
        anchor = (0.5, -0.2)
    else:
        ax.set_ylabel("Liquid hydrocarbon LHV\n[MW]", fontsize=24, weight="bold")
        ax.set_xlabel("Shale region", fontsize=24, weight="bold")
        ax.set_ylim(0, 950)
        ax.tick_params(axis="x", rotation=45)
        ax.bar_label(bars, label_type="edge", fontsize=16, weight="bold", fmt="%.0f", rotation=90)
        anchor = (0.5, -0.3)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(loc="upper center", bbox_to_anchor=anchor, ncol=4, fontsize=20)
    _save(fig, path)


def _plot_remaining(
    output_dir: Path, rok: dict[str, Any], regional: dict[str, Any]
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 20})
    generated: dict[str, Any] = OrderedDict()

    _stacked_emissions(
        ["Literature", "Optimal\nwith M2", "Optimal\nwith M3", "Optimal\nwith M4", "Optimal\nwith M5"],
        [9.26] + [rok[name]["upstream_emissions_g_co2e_per_mj"] for name in ROK_MODELS],
        [14.2] + [rok[name]["process_emissions_g_co2e_per_mj"] for name in ROK_MODELS],
        output_dir / OUTPUT_NAMES["F3"],
    )
    generated["F3"] = OUTPUT_NAMES["F3"]

    tax_rows = _read_csv(REPO_ROOT / "optimal_data_wrt_c_tax_rates.csv")
    fig, ax = plt.subplots(figsize=(8, 6))
    x = [float(row["C-tax-rate"]) for row in tax_rows]
    # The notebook uses 0.01 USD/tonne as its smallest positive log-scale case;
    # its zero-tax optimum is numerically redundant and cannot appear at log(0).
    positive = [(tax, float(row["MSP"])) for tax, row in zip(x, tax_rows) if tax > 0]
    ax.semilogx([p[0] for p in positive], [p[1] for p in positive], marker="o", markersize=10, linewidth=3)
    ax.set_xlabel("CO$_2$ emissions tax rate\n[USD/tonne CO$_2$e]", fontsize=24, weight="bold")
    ax.set_ylabel("Minimum Selling Price\n[USD/MJ fuel]", fontsize=24, weight="bold")
    ax.grid(which="major")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(axis="both", labelsize=16, length=6, which="both")
    _save(fig, output_dir / OUTPUT_NAMES["F5"])
    generated["F5"] = OUTPUT_NAMES["F5"]

    _stacked_emissions(
        list(REGIONS),
        [regional[r]["upstream_emissions_g_co2e_per_mj"] for r in REGIONS],
        [regional[r]["process_emissions_g_co2e_per_mj"] for r in REGIONS],
        output_dir / OUTPUT_NAMES["F6"],
    )
    generated["F6"] = OUTPUT_NAMES["F6"]
    _component_bars(regional, REGIONS, "liquid_component_lhv_mw", False, output_dir / OUTPUT_NAMES["F7"])
    generated["F7"] = OUTPUT_NAMES["F7"]

    regional_rows = {row["Region"]: row for row in _read_csv(REPO_ROOT / "optimal_data_wrt_region.csv")}
    ng_fraction = _read_csv(REPO_ROOT / "data" / "NGL_fraction.csv")[0]
    methane_percent = {r: 100.0 * (1.0 - float(ng_fraction[r])) for r in REGIONS}
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter([methane_percent[r] for r in REGIONS], [float(regional_rows[r]["MSP"]) for r in REGIONS], s=60)
    for index, region in enumerate(REGIONS):
        offset = (-50, 10) if region == "EF-Basin" else (-50, -15) if region == "Bakken" else ((-40, 0) if index % 2 == 0 else (5, 0))
        ax.annotate(region, (methane_percent[region], float(regional_rows[region]["MSP"])), xytext=offset, textcoords="offset points", fontsize=14)
    ax.set_ylabel("Minimum Selling Price\n[USD/MJ fuel]", fontsize=24, weight="bold")
    ax.set_xlabel("Mole-percent of CH$_4$ in shale [%]", fontsize=24, weight="bold")
    ax.set_xlim(40, 100)
    ax.grid()
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(axis="both", labelsize=16)
    _save(fig, output_dir / OUTPUT_NAMES["F8"])
    generated["F8"] = OUTPUT_NAMES["F8"]

    _component_bars(rok, ROK_MODELS, "liquid_component_flow_mol_per_s", True, output_dir / OUTPUT_NAMES["FS1"])
    generated["FS1"] = OUTPUT_NAMES["FS1"]
    fig, ax = plt.subplots(figsize=(8, 6))
    values = [rok[name]["c4_plus_olefin_mole_percent"] for name in ROK_MODELS]
    bars = ax.bar(ROK_MODELS, values)
    ax.set_ylabel("C$_{4+}$ olefin outlet\nmole-percent [%]", fontsize=24, weight="bold")
    ax.set_xlabel("ROK model", fontsize=24, weight="bold")
    ax.set_ylim(0, 100)
    ax.bar_label(bars, label_type="edge", weight="bold", fontsize=16, fmt="%.2f")
    ax.tick_params(axis="y", labelsize=20)
    _save(fig, output_dir / OUTPUT_NAMES["FS2"])
    generated["FS2"] = OUTPUT_NAMES["FS2"]

    shale_rows = _read_csv(REPO_ROOT / "data" / "Shale_composition.csv")
    feed = OrderedDict((r, OrderedDict()) for r in REGIONS)
    for row in shale_rows:
        for region in REGIONS:
            feed[region][row["Species"]] = 100.0 * float(row[region])
    fig, ax = plt.subplots(figsize=(8, 6))
    left = np.zeros(len(REGIONS))
    plotted_components = [
        component
        for component in feed["EF-1"]
        if any(feed[region][component] != 0.0 for region in REGIONS)
    ]
    for component in plotted_components:
        vals = [feed[r][component] for r in REGIONS]
        ax.barh(REGIONS, vals, left=left, label=component, color=colors_dict[component])
        left += vals
    ax.set_xlabel("Feed mole percent [%]", fontsize=24, weight="bold")
    ax.set_ylabel("Region", fontsize=24, weight="bold")
    ax.set_xlim(0, 100)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=20)
    _save(fig, output_dir / OUTPUT_NAMES["FS5"])
    generated["FS5"] = OUTPUT_NAMES["FS5"]
    return generated


def _reference_differences(rok: dict[str, Any], regional: dict[str, Any]) -> dict[str, Any]:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    f3_ref = {row["case"]: row for row in reference["figure_3"]["series"]}
    f6_ref = {row["region"]: row for row in reference["figure_6"]["series"]}
    return {
        "F3": {
            name: {
                "upstream_difference": rok[name]["upstream_emissions_g_co2e_per_mj"] - f3_ref[name]["upstream"],
                "process_difference": rok[name]["process_emissions_g_co2e_per_mj"] - f3_ref[name]["process"],
                "known_issue": "EMISSIONS-NORMALIZATION-001",
            }
            for name in ROK_MODELS
        },
        "F6": {
            region: {
                "upstream_difference": regional[region]["upstream_emissions_g_co2e_per_mj"] - f6_ref[region]["upstream"],
                "process_difference": regional[region]["process_emissions_g_co2e_per_mj"] - f6_ref[region]["process"],
                "known_issue": "EMISSIONS-NORMALIZATION-001",
            }
            for region in REGIONS
        },
        "F7_total_difference_mw": {
            region: regional[region]["total_liquid_lhv_mw"] - reference["figure_7"]["total_labels"][region]
            for region in REGIONS
        },
        "FS1_total_difference_mol_per_s": {
            name: rok[name]["total_liquid_flow_mol_per_s"] - reference["figure_s1"]["total_labels"][name]
            for name in ROK_MODELS
        },
        "FS2_difference_mole_percent": {
            name: rok[name]["c4_plus_olefin_mole_percent"] - reference["figure_s2"]["values"][name]
            for name in ROK_MODELS
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--record", required=True, type=Path)
    args = parser.parse_args()
    if REPO_ROOT == args.output_dir.resolve() or REPO_ROOT in args.output_dir.resolve().parents:
        raise ValueError("Output directory must be outside the repository")

    rok, regional, hashes = collect_checkpoint_data()
    generated = _plot_remaining(args.output_dir, rok, regional)
    report = {
        "schema_version": 1,
        "experiment_id": "A1-ARCHIVED-NONCOMPOSITE-FIGURES-001",
        "status": "complete_with_classified_published_snapshot_differences",
        "coverage": list(OUTPUT_NAMES),
        "source_hashes": hashes,
        "csv_source_hashes": {
            str(path.relative_to(REPO_ROOT)): _sha256(path)
            for path in (
                REPO_ROOT / "optimal_data_wrt_c_tax_rates.csv",
                REPO_ROOT / "optimal_data_wrt_region.csv",
                REPO_ROOT / "data" / "NGL_fraction.csv",
                REPO_ROOT / "data" / "Shale_composition.csv",
            )
        },
        "generated_pdfs": generated,
        "rok_models": rok,
        "regions": regional,
        "published_reference_differences": _reference_differences(rok, regional),
        "output_policy": "Generated PDFs and full numeric record are outside the repository; archived PDFs were not modified.",
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    args.record.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "coverage": report["coverage"], "output_dir": str(args.output_dir), "record": str(args.record)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
