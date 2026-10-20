"""Verify a sealed Sections 4.4--4.6 direct-statistics evidence bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


FORBIDDEN_OUTPUTS = frozenset(
    {
        "oof_predictions.csv",
        "oof_prediction_partitions.csv",
        "fold_metrics.csv",
        "fold_coefficients.csv",
        "model_summary.csv",
        "multiplicity_decisions.csv",
        "analysis_control_manifest.csv",
        "matched_ess_diagnostics.csv",
        "weight_permutation_diagnostics.csv",
    }
)
REQUIRED = (
    "gate.json",
    "RUN_REPORT.md",
    "RESULTS_INDEX.md",
    "evidence_manifest.json",
    "sha256_manifest.csv",
    "cpp_feature_map_manifest.csv",
    "direct_effects.csv",
    "single_map_advantages.csv",
    "weighting_effects.csv",
    "bootstrap_summary.csv",
    "support_summary.csv",
    "common_support_manifest.csv",
    "weighting_mechanism_diagnostics.csv",
    "report_artifact_contract.csv",
    "time_projection.json",
)
CUMULANT_CENTERING_VERSION = "order_ge_2.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _safe_member(root: Path, member: str) -> Path:
    candidate = (root / member).resolve()
    if candidate != root and root not in candidate.parents:
        raise ValueError(f"Manifest path escapes the run directory: {member}")
    return candidate


def _verify_sha256_manifest(root: Path) -> tuple[bool, list[str]]:
    manifest_path = root / "sha256_manifest.csv"
    failures: list[str] = []
    with manifest_path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return False, ["sha256 manifest has no archived files"]
    for row in rows:
        relative = str(row.get("path", ""))
        expected = str(row.get("sha256", ""))
        target = _safe_member(root, relative)
        if not target.is_file():
            failures.append(f"missing:{relative}")
        elif _sha256(target) != expected:
            failures.append(f"hash_mismatch:{relative}")
    return not failures, failures


def _verify_current_source_hashes(manifest: dict[str, Any]) -> list[str]:
    environment = manifest.get("environment", {})
    sources = environment.get("source_files", []) if isinstance(environment, dict) else []
    failures: list[str] = []
    for item in sources:
        source = Path(str(item.get("path", "")))
        expected = str(item.get("sha256", ""))
        if not source.is_file():
            failures.append(f"source_missing:{source}")
        elif _sha256(source) != expected:
            failures.append(f"source_hash_mismatch:{source}")
    return failures


def _metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def _verify_cumulant_centering(root: Path, expected_maps: int) -> list[str]:
    manifest_path = root / "cpp_feature_map_manifest.csv"
    with manifest_path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    failures: list[str] = []
    required = {
        "case_id",
        "order",
        "algorithm",
        "metadata_path",
        "metadata_sha256",
        "cumulant_centering",
        "cumulant_centering_version",
    }
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        return [f"cpp_manifest_missing_columns:{','.join(sorted(missing))}"]
    if len(rows) != expected_maps:
        failures.append(f"cpp_manifest_rows:{len(rows)}")
    for row in rows:
        case_id = str(row["case_id"])
        try:
            order = int(row["order"])
        except ValueError:
            failures.append(f"cpp_manifest_order:{case_id}")
            continue
        expected_centering = "none" if order == 1 else (
            "weighted_window_mean"
            if str(row["algorithm"]) == "Weighted moment"
            else "window_mean"
        )
        if str(row["cumulant_centering"]) != expected_centering:
            failures.append(f"cpp_manifest_centering:{case_id}")
        if str(row["cumulant_centering_version"]) != CUMULANT_CENTERING_VERSION:
            failures.append(f"cpp_manifest_centering_version:{case_id}")
        metadata_path = Path(str(row["metadata_path"])).resolve()
        if root != metadata_path and root not in metadata_path.parents:
            failures.append(f"cpp_metadata_outside_bundle:{case_id}")
            continue
        if not metadata_path.is_file():
            failures.append(f"cpp_metadata_missing:{case_id}")
            continue
        if _sha256(metadata_path) != str(row["metadata_sha256"]):
            failures.append(f"cpp_metadata_hash:{case_id}")
            continue
        metadata = _metadata(metadata_path)
        if metadata.get("cumulant_centering") != expected_centering:
            failures.append(f"cpp_metadata_centering:{case_id}")
        if metadata.get("cumulant_centering_version") != CUMULANT_CENTERING_VERSION:
            failures.append(f"cpp_metadata_centering_version:{case_id}")
        if metadata.get("center_cumulant_input_for_order_ge_2") != "true":
            failures.append(f"cpp_metadata_centering_request:{case_id}")
    return failures


def verify_bundle(
    run_dir: Path,
    *,
    mode: str,
    require_current_sources: bool,
    require_time_budget: bool,
) -> tuple[bool, list[str]]:
    root = Path(run_dir).resolve()
    failures = [f"missing:{name}" for name in REQUIRED if not (root / name).is_file()]
    if failures:
        return False, failures
    gate = _read_json(root / "gate.json")
    evidence = _read_json(root / "evidence_manifest.json")
    expected_status = "PASS" if mode == "formal" else "ANALYSIS_SMOKE_PASS"
    if gate.get("status") != expected_status:
        failures.append(f"gate_status:{gate.get('status')}")
    if gate.get("mode") != mode:
        failures.append(f"gate_mode:{gate.get('mode')}")
    expected_maps = 102 if mode == "formal" else 34
    if int(gate.get("map_counts", {}).get("map_count", -1)) != expected_maps:
        failures.append(f"map_count:{gate.get('map_counts', {})}")
    failures.extend(_verify_cumulant_centering(root, expected_maps))
    failed_checks = [str(item.get("name", "unnamed")) for item in gate.get("checks", []) if item.get("passed") is not True]
    failures.extend(f"failed_check:{name}" for name in failed_checks)
    manifest_ok, manifest_failures = _verify_sha256_manifest(root)
    if not manifest_ok:
        failures.extend(manifest_failures)
    forbidden = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if relative.parts[:1] == ("section_4_4",) and path.name.lower() == "multiplicity_decisions.csv":
            continue
        name = path.name.lower()
        if name in FORBIDDEN_OUTPUTS or (name.startswith("learned_threshold") and name.endswith(".csv")):
            forbidden.append(relative.as_posix())
    failures.extend(f"forbidden:{path}" for path in forbidden)
    projection = _read_json(root / "time_projection.json")
    formal = projection.get("formal_projection", {})
    if require_time_budget and formal.get("time_budget_pass") is not True:
        failures.append("formal_time_budget_projection:fail")
    if require_current_sources:
        failures.extend(_verify_current_source_hashes(evidence))
    return not failures, failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--mode", choices=("analysis-smoke", "formal"), required=True)
    parser.add_argument("--require-current-sources", action="store_true")
    parser.add_argument("--require-time-budget", action="store_true")
    arguments = parser.parse_args(argv)
    passed, failures = verify_bundle(
        arguments.run_dir,
        mode=arguments.mode,
        require_current_sources=arguments.require_current_sources,
        require_time_budget=arguments.require_time_budget,
    )
    print("bundle_verification=" + ("PASS" if passed else "FAIL"))
    print("failures=" + (",".join(failures) if failures else "none"))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
