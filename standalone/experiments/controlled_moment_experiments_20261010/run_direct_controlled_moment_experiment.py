"""Run classifier-free matching-moment controls for the direct protocol."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from controlled_distributions import controlled_experiments
from direct_evaluation import (
    DirectControlledConfig,
    direct_control_gates,
    evaluate_direct_controls,
    finalize_direct_control_gates,
    summarize_direct_controls,
)


EXPERIMENT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_DIRECTORY = EXPERIMENT_DIRECTORY.parents[1]
DEFAULT_OUTPUT_ROOT = REPOSITORY_DIRECTORY / "output" / "direct_controlled_moment_experiments_20261011"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_config(path: Path) -> dict[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Direct controlled configuration must contain a JSON object")
    return raw


def _quick_config(config: DirectControlledConfig) -> DirectControlledConfig:
    return DirectControlledConfig(
        random_seeds=config.random_seeds[:2],
        window_sample_counts=config.window_sample_counts,
        windows_per_class=min(config.windows_per_class, 100),
        high_threshold=config.high_threshold,
        variance_epsilon=config.variance_epsilon,
        bootstrap_replicates=min(config.bootstrap_replicates, 200),
    )


def _write_report(path: Path, summary, gates: dict[str, object], run_kind: str) -> None:
    lines = [
        "# Direct Matching-Moment Control Report",
        "",
        f"Run type: `{run_kind}`.",
        "",
        "The control compares finite-window cumulant estimators directly. It does not fit a classifier, split training and test samples, or compute predictive metrics.",
        "",
        "## Target-order estimates",
        "",
        "| Pair | Samples | Order | Analytic difference | Estimated difference | 95% CI | Detected |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    targets = summary.loc[summary["order"] == summary["target_order"]]
    for row in targets.itertuples(index=False):
        lines.append(
            f"| {row.experiment} | {row.window_samples} | {row.order} | {row.analytic_delta_a_minus_b:.8g} | "
            f"{row.estimate_delta_mean:.8g} | [{row.estimate_delta_ci95_low:.8g}, {row.estimate_delta_ci95_high:.8g}] | "
            f"{row.target_order_detected} |"
        )
    lines.extend(["", "## Gates", "", "| Gate | Value |", "| --- | --- |"])
    for key, value in gates.items():
        lines.append(f"| {key} | {value} |")
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `population_moments.csv`: analytic cumulants for the paired distributions.",
            "- `direct_per_seed_estimates.csv`: finite-window direct cumulant estimates by seed.",
            "- `direct_summary.csv`: bias, RMSE, group differences, and confidence intervals.",
            "- `gate.json`: execution and target-order detection gates.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(config_path: Path, output_root: Path, *, quick: bool = False, output_directory: Path | None = None) -> Path:
    config_path = config_path.resolve()
    raw_config = _load_config(config_path)
    config = DirectControlledConfig.from_mapping(raw_config)
    if quick:
        config = _quick_config(config)
    run_kind = "smoke" if quick else "formal"
    run_dir = output_directory.resolve() if output_directory is not None else output_root.resolve() / datetime.now(UTC).strftime("run_%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=False)
    experiments = controlled_experiments()
    population_rows = [row for experiment in experiments for row in experiment.moment_rows()]
    results = evaluate_direct_controls(experiments, config)
    summary = summarize_direct_controls(results, config)
    gates = finalize_direct_control_gates(direct_control_gates(experiments, results, summary))
    gates["run_kind"] = run_kind
    gates["component_execution"] = "COMPLETE" if bool(gates["overall_pass"]) else "FAILED"
    gates["status"] = "VALIDATION_PASS" if bool(gates["overall_pass"]) else "FAIL"
    population_path = run_dir / "population_moments.csv"
    results_path = run_dir / "direct_per_seed_estimates.csv"
    summary_path = run_dir / "direct_summary.csv"
    population_path.write_text(__import__("pandas").DataFrame(population_rows).to_csv(index=False), encoding="utf-8")
    results.to_csv(results_path, index=False)
    summary.to_csv(summary_path, index=False)
    (run_dir / "run_config.json").write_text(
        json.dumps({"run_kind": run_kind, **config.as_jsonable()}, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (run_dir / "gate.json").write_text(json.dumps(gates, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(run_dir / "RUN_REPORT.md", summary, gates, run_kind)
    source_files = [
        EXPERIMENT_DIRECTORY / "controlled_distributions.py",
        EXPERIMENT_DIRECTORY / "direct_evaluation.py",
        Path(__file__),
        config_path,
    ]
    output_files = sorted(path for path in run_dir.iterdir() if path.is_file())
    manifest = {
        "schema_version": "direct-controlled-moment-evidence.v1",
        "generated_utc": datetime.now(UTC).isoformat(),
        "run_kind": run_kind,
        "config": config.as_jsonable(),
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "gates": gates,
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
        "source_files": [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in source_files],
        "output_files": [{"path": path.name, "sha256": _sha256(path)} for path in output_files],
    }
    (run_dir / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return run_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run direct matching-moment cumulant controls without classifiers.")
    parser.add_argument("--config", type=Path, default=EXPERIMENT_DIRECTORY / "experiment_config.json")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--quick", action="store_true")
    arguments = parser.parse_args(argv)
    result = run(arguments.config, arguments.output_root, quick=arguments.quick, output_directory=arguments.output_directory)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
