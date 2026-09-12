#!/usr/bin/env python3
"""Export the M5/Bakken NL file and its actual Pyomo symbol ordering."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
import weakref
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import idaes
import pyomo
from pyomo.opt import ProblemFormat

from reproducibility.run_m5_bakken_tax_series import (
    FileDeterminism,
    _build_preoptimization_model,
    _load_column_order,
    _write_report,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


SYMBOL = re.compile(r"^([vco])(\d+)$")


def _resolve_symbol_object(item: Any) -> Any:
    return item() if isinstance(item, weakref.ReferenceType) else item


def _ordered_symbol_names(symbol_map: Any) -> dict[str, list[str]]:
    grouped: dict[str, list[tuple[int, str]]] = {
        "variables": [],
        "constraints": [],
        "objectives": [],
    }
    categories = {"v": "variables", "c": "constraints", "o": "objectives"}
    for symbol, reference in symbol_map.bySymbol.items():
        match = SYMBOL.match(symbol)
        if match is None:
            continue
        component = _resolve_symbol_object(reference)
        if component is None:
            raise RuntimeError(f"Symbol {symbol} points to an expired weak reference.")
        grouped[categories[match.group(1)]].append(
            (int(match.group(2)), component.name)
        )
    return {
        category: [name for _, name in sorted(entries)]
        for category, entries in grouped.items()
    }


def _digest(names: list[str]) -> str:
    return hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()


def _inspect_nl(path: Path, header_lines: int = 10) -> dict[str, Any]:
    digest = hashlib.sha256()
    header: list[str] = []
    size = 0
    with path.open("rb") as stream:
        for _ in range(header_lines):
            line = stream.readline()
            if not line:
                break
            digest.update(line)
            size += len(line)
            header.append(line.decode(errors="replace").rstrip("\r\n"))
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
    return {
        "path": str(path),
        "bytes": size,
        "sha256": digest.hexdigest(),
        "header": header,
    }


def export_symbol_map(
    nl_output: Path,
    determinism: str,
    column_order_path: Path | None = None,
    inline_defined_variables: bool = False,
) -> dict[str, Any]:
    model = _build_preoptimization_model(5, "Bakken", 0.0)
    unfix_DOFs_pre_optimization(model)
    option = {
        "ordered": FileDeterminism.ORDERED,
        "sort-indices": FileDeterminism.SORT_INDICES,
        "sort-symbols": FileDeterminism.SORT_SYMBOLS,
    }[determinism]
    nl_output.parent.mkdir(parents=True, exist_ok=True)
    io_options = {
        "file_determinism": option,
        # The in-memory SymbolMap retains component objects without embedding
        # their very long names on every NL expression line.
        "symbolic_solver_labels": False,
    }
    if inline_defined_variables:
        io_options["export_defined_variables"] = False
    requested_column_order_sha256 = None
    if column_order_path is not None:
        column_order, requested_column_order_sha256 = _load_column_order(
            model, column_order_path
        )
        io_options["column_order"] = column_order
    filename, symbol_map_id = model.write(
        str(nl_output),
        format=ProblemFormat.nl,
        io_options=io_options,
    )
    symbol_map = model.solutions.symbol_map[symbol_map_id]
    order = _ordered_symbol_names(symbol_map)
    nl_path = Path(filename)
    return {
        "schema_version": 1,
        "experiment_id": "B-M5-BAKKEN-NL-SYMBOL-MAP-001",
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "idaes": idaes.__version__,
            "pyomo": pyomo.__version__,
        },
        "case": {
            "model_code": 5,
            "region": "Bakken",
            "co2_tax_usd_per_kg": 0.0,
            "free_design_variables": 8,
            "file_determinism": determinism,
            "symbolic_solver_labels": False,
            "inline_defined_variables": inline_defined_variables,
            "column_order_from_symbol_map": (
                str(column_order_path.resolve())
                if column_order_path is not None
                else None
            ),
            "requested_column_order_sha256": requested_column_order_sha256,
        },
        "nl_file": _inspect_nl(nl_path),
        "ordering": {
            category: {
                "count": len(names),
                "sha256": _digest(names),
                "names": names,
            }
            for category, names in order.items()
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nl-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--file-determinism",
        choices=("ordered", "sort-indices", "sort-symbols"),
        default="sort-symbols",
    )
    parser.add_argument("--column-order-from-symbol-map", type=Path)
    parser.add_argument(
        "--inline-defined-variables",
        action="store_true",
        help=(
            "Set export_defined_variables=false so modern Pyomo inlines "
            "named Expression objects instead of emitting NL common expressions."
        ),
    )
    arguments = parser.parse_args()
    report = export_symbol_map(
        arguments.nl_output,
        arguments.file_determinism,
        arguments.column_order_from_symbol_map,
        arguments.inline_defined_variables,
    )
    _write_report(report, arguments.output)
    print(
        json.dumps(
            {
                "environment": report["environment"],
                "case": report["case"],
                "nl_file": report["nl_file"],
                "ordering": {
                    category: {
                        "count": details["count"],
                        "sha256": details["sha256"],
                    }
                    for category, details in report["ordering"].items()
                },
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
