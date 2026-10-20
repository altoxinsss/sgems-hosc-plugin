from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from run_direct_controlled_moment_experiment import run


def test_direct_control_quick_run_writes_only_direct_cumulant_evaluation_outputs(tmp_path: Path) -> None:
    experiment_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "direct_control_smoke"
    result = run(
        experiment_root / "experiment_config.json",
        tmp_path,
        quick=True,
        output_directory=output,
    )
    assert result == output
    summary = pd.read_csv(output / "direct_summary.csv")
    estimates = pd.read_csv(output / "direct_per_seed_estimates.csv")
    prohibited_columns = {"classifier", "model", "fold", "oof", "pr_auc", "f1", "threshold"}
    assert prohibited_columns.isdisjoint(summary.columns.str.lower())
    assert prohibited_columns.isdisjoint(estimates.columns.str.lower())
    gate = json.loads((output / "gate.json").read_text(encoding="utf-8"))
    assert gate["run_kind"] == "smoke"
    assert set(summary["statistic"]) == {"kappa1", "kappa2", "kappa3", "kappa4", "kappa5", "kappa6"}
