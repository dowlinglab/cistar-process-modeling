#!/usr/bin/env python3
"""Export a complete named-variable state from an M5/Bakken checkpoint."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from idaes.core.util import model_serializer as ms

from reproducibility.run_m5_bakken_tax_series import (
    _build_preoptimization_model,
    _capture_named_state,
    _checkpoint,
)
from src.unit_initialization import unfix_DOFs_pre_optimization


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archived-optimal-tax", type=float)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    model = _build_preoptimization_model()
    tax_rate = 0.0
    stage = "historical-costed-checkpoint"
    if arguments.archived_optimal_tax is not None:
        tax_rate = arguments.archived_optimal_tax
        checkpoint = _checkpoint(
            "CISTAR_optimal_solution_Bakken_C_tax_{}_M5_purge_0.01_"
            "sequential_solve.json.gz".format(tax_rate)
        )
        if not checkpoint.exists():
            raise FileNotFoundError(f"Archived starting point not found: {checkpoint}")
        model.fs.c_tax_rate = tax_rate
        unfix_DOFs_pre_optimization(model)
        ms.from_json(model, fname=str(checkpoint))
        # Serializer checkpoints preserve fixed statuses. Match the solver-start
        # state used by the tax runner after loading an archived optimum.
        unfix_DOFs_pre_optimization(model)
        stage = f"archived-optimal-tax-{tax_rate}"
    summary = _capture_named_state(
        model,
        arguments.output,
        stage=stage,
        tax_rate=tax_rate,
    )
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
