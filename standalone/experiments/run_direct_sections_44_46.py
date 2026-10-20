"""Formal Sections 4.4--4.6 direct-statistics experiment entrypoint.

This entrypoint is separate from the historical 486-map model protocol.  It
never trains a classifier or writes fold, threshold, coefficient, or OOF
artifacts.  Section 4.4 remains a template-based spatial-cumulant component;
Sections 4.5 and 4.6 compare local marginal maps directly to channel facies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


EXPERIMENTS_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = EXPERIMENTS_ROOT.parent
R15_ROOT = EXPERIMENTS_ROOT / "r1_5_incremental_value_20260924"
PREPROCESS_ROOT = EXPERIMENTS_ROOT / "stanford_vie_preprocessing"
TEMPLATE_ROOT = EXPERIMENTS_ROOT / "template_directional_20261010"
CONTROLLED_ROOT = EXPERIMENTS_ROOT / "controlled_moment_experiments_20261010"
CONFIG_ROOT = EXPERIMENTS_ROOT / "configs"
for directory in (R15_ROOT, CONTROLLED_ROOT):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from cli_runner import (
    CUMULANT_CENTERING_VERSION,
    CppFeatureMapArtifact,
    CppFeatureMapSpec,
    build_direct_cpp_feature_map_specs,
    run_cpp_feature_matrix,
    sha256,
    validate_direct_cpp_feature_map_protocol,
)
from cpp_python_parity_gate import run_cpp_python_fixture_matrix_gate, run_cpp_python_parity_gate, run_uniform_three_route_audit
from direct_statistical_evaluation import (
    DEFAULT_BASELINES,
    DEFAULT_HIGH_FEATURES,
    DirectCase,
    apply_common_support,
    annotate_signed_direction_consistency,
    build_analysis_support,
    compute_direct_effects,
    compute_parameter_sensitivity_effects,
    compute_single_map_advantages,
    compute_weighting_effects,
    derive_standardized_cumulants,
    joint_spatial_block_bootstrap,
)
from evidence import collect_environment, seal_evidence_bundle, write_json
from grid_io import GridData, read_grid
from mask_protocol import write_mask_bundle
from report_direct_results import generate_direct_report
from weighting import (
    WeightConfig,
    compute_weighted_statistics,
    coordinate_gaussian_window_weights,
    generate_weights,
)
from window_support import extract_complete_windows


DIRECT_REQUIRED_ARTIFACTS = (
    "cpp_feature_map_manifest.csv",
    "direct_feature_summary.csv",
    "direct_effects.csv",
    "single_map_advantages.csv",
    "weighting_effects.csv",
    "bootstrap_summary.csv",
    "support_summary.csv",
    "direct_analysis_values.csv",
    "common_support_manifest.csv",
    "parameter_sensitivity_contrasts.csv",
    "lambda4_independent_audit.csv",
    "weighting_mechanism_diagnostics.csv",
    "time_projection.json",
    "report_artifact_contract.csv",
    "RUN_REPORT.md",
    "RESULTS_INDEX.md",
    "gate.json",
    "evidence_manifest.json",
    "sha256_manifest.csv",
)
# The direct protocol may retain ordinary terms such as ``threshold`` in
# configuration provenance and component reports.  Only the historical
# model-protocol artifacts themselves are prohibited from this run bundle.
FORBIDDEN_OUTPUT_ARTIFACTS = frozenset(
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
FORBIDDEN_OUTPUT_PREFIXES = ("learned_threshold",)
FEATURES = (*DEFAULT_BASELINES, *DEFAULT_HIGH_FEATURES)
ANALYSIS_SMOKE_BOOTSTRAP_REPLICATES = 100


def _timestamp() -> str:
    return datetime.now(UTC).strftime("run_%Y%m%d_%H%M%S")


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Configuration must contain an object: {path}")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _resolve(config_path: Path, value: str | Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (config_path.parent / path).resolve()


def _combined_config_hash(paths: Mapping[str, Path]) -> str:
    digest = hashlib.sha256()
    for name, path in sorted(paths.items()):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _run_command(command: Sequence[str], log_path: Path, *, cwd: Path, env: Mapping[str, str] | None = None) -> dict[str, Any]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    # Build tools on this Windows host can emit text in a console encoding that
    # differs from Python's default.  A diagnostic log must not turn that
    # nonessential decoding difference into a failed experiment run.
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
        env=env,
    )
    log_path.write_text(
        "command=" + json.dumps(list(command)) + "\n\nstdout\n" + completed.stdout + "\n\nstderr\n" + completed.stderr,
        encoding="utf-8",
    )
    return {"command": list(command), "exit_code": int(completed.returncode), "log_path": str(log_path.resolve())}


def _require_command(record: Mapping[str, Any], name: str) -> dict[str, Any]:
    if int(record["exit_code"]) != 0:
        raise RuntimeError(f"{name} failed; see {record['log_path']}")
    return {"name": name, "passed": True, "details": str(record["log_path"])}


def _create_run(output_root: Path) -> Path:
    run_dir = output_root.resolve() / _timestamp()
    suffix = 1
    while run_dir.exists():
        run_dir = output_root.resolve() / f"{_timestamp()}_{suffix:02d}"
        suffix += 1
    for relative in (
        "configs", "logs", "preprocessing", "section_4_4", "section_4_5/cpp_maps", "section_4_5/derived_maps",
        "section_4_6/cpp_maps", "section_4_6/sensitivity", "maps/cpp", "tables", "figures", "support_masks", "components",
    ):
        (run_dir / relative).mkdir(parents=True, exist_ok=True)
    return run_dir


def _copy_frozen_configs(run_dir: Path, sources: Mapping[str, Path]) -> dict[str, Path]:
    frozen: dict[str, Path] = {}
    for name, source in sources.items():
        # Several component configurations are named experiment_config.json.
        # The logical source key is part of the frozen filename so every input
        # remains recoverable in the sealed evidence bundle.
        target = run_dir / "configs" / f"{name}{source.suffix}"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        frozen[name] = target
    return frozen


def _build_release_and_tests(python: Path, common: Mapping[str, Any], run_dir: Path, *, common_path: Path | None = None) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    common_path = common_path or CONFIG_ROOT / "direct_common_protocol.json"
    build_dir = _resolve(common_path, str(common["cmake_build_dir"]))
    ctest_value = str(common.get("ctest_executable", "ctest"))
    cmake_value = str(common.get("cmake_executable", "cmake"))
    ctest = Path(shutil.which(ctest_value) or _resolve(common_path, ctest_value))
    cmake = Path(shutil.which(cmake_value) or _resolve(common_path, cmake_value))
    if not build_dir.is_dir() or not ctest.is_file() or not cmake.is_file():
        raise FileNotFoundError("Direct protocol Release build or CTest executable is unavailable")
    checks.append(_require_command(
        _run_command([str(cmake), "--build", str(build_dir), "--config", "Release", "--parallel"], run_dir / "logs" / "cmake_release_build.log.txt", cwd=REPOSITORY_ROOT),
        "release_build",
    ))
    test_paths = [_resolve(common_path, value) for value in common["python_test_paths"]]
    test_environment = dict(os.environ)
    test_environment["WITHOUT_SGEMS_BUILD_DIR"] = str(_resolve(common_path, str(common["movewindow_cli"])).parent)
    checks.append(_require_command(
        _run_command([str(python), "-m", "pytest", *[str(path) for path in test_paths]], run_dir / "logs" / "python_tests.log.txt", cwd=EXPERIMENTS_ROOT, env=test_environment),
        "python_regression_tests",
    ))
    checks.append(_require_command(
        _run_command([str(ctest), "--test-dir", str(build_dir), "-C", "Release", "--output-on-failure"], run_dir / "logs" / "ctest.log.txt", cwd=REPOSITORY_ROOT),
        "native_release_ctest",
    ))
    return checks


def _preprocessing_runtime_config(source: Path, output_root: Path) -> dict[str, Any]:
    """Freeze preprocessing inputs before moving its configuration into a run."""

    runtime = _load_json(source)
    for key in ("porosity_path", "facies_path", "poro_inc_path", "co_location_fixture_path"):
        if key in runtime:
            runtime[key] = str(_resolve(source, str(runtime[key])))
    runtime["output_root"] = str(output_root.resolve())
    return runtime


def _run_preprocessing(python: Path, common_source: Path, common: Mapping[str, Any], run_dir: Path) -> dict[str, Any]:
    source = _resolve(common_source, str(common["preprocessing_config"]))
    runtime = _preprocessing_runtime_config(source, run_dir / "preprocessing")
    config_path = run_dir / "configs" / "preprocessing_runtime.json"
    _write_json(config_path, runtime)
    record = _run_command(
        [str(python), str(PREPROCESS_ROOT / "prepare_stanford_vie.py"), "--config", str(config_path)],
        run_dir / "logs" / "preprocessing.log.txt",
        cwd=PREPROCESS_ROOT,
    )
    _require_command(record, "stanford_preprocessing")
    manifest_path = run_dir / "preprocessing" / "preprocessing_manifest.json"
    manifest = _load_json(manifest_path)
    if float(manifest["poro_inc_max_absolute_error"]) > float(runtime["max_poro_inc_absolute_error"]):
        raise RuntimeError("Preprocessing poro.inc co-location audit did not pass")
    if len(manifest.get("subvolumes", [])) != 3:
        raise RuntimeError("Direct protocol preprocessing did not create three fixed subvolumes")
    return {"name": "stanford_preprocessing", "passed": True, "details": str(manifest_path.resolve()), "manifest": manifest}


def _load_direct_grids(common_path: Path, common: Mapping[str, Any], preprocessing_root: Path) -> tuple[dict[str, GridData], dict[str, dict[str, GridData]], dict[str, int]]:
    porosity: dict[str, GridData] = {}
    references: dict[str, dict[str, GridData]] = {}
    representative_z: dict[str, int] = {}
    for item in common["datasets"]:
        dataset_id = str(item["id"])
        source_path = _resolve(common_path, str(item["porosity"]))
        name = source_path.name
        local_path = preprocessing_root / dataset_id / name
        if not local_path.is_file():
            raise FileNotFoundError(f"Run-local preprocessing map is missing: {local_path}")
        porosity[dataset_id] = read_grid(local_path, "porosity")
        references[dataset_id] = {}
        for key, property_name in (("label", "channel_strict_majority"), ("label_valid", "channel_label_valid"), ("occupancy", "channel_occupancy")):
            configured = _resolve(common_path, str(item[key]))
            candidate = preprocessing_root / dataset_id / configured.name
            references[dataset_id][key] = read_grid(candidate, property_name)
        representative_z[dataset_id] = int(item["representative_z"])
    return porosity, references, representative_z


def _runtime_template_config(common_path: Path, common: Mapping[str, Any], frozen_section44: Path, preprocessing_root: Path, run_dir: Path, entrypoint_id: str) -> Path:
    source = _resolve(common_path, str(common["template_component_config"]))
    config = _load_json(source)
    for key in ("template_cli", "native_test_workdir"):
        config[key] = str(_resolve(source, str(config[key])))
    config["output_root"] = str((run_dir / "section_4_4").resolve())
    config["preprocessing_manifest"] = str((preprocessing_root / "preprocessing_manifest.json").resolve())
    config["formal_section44_protocol"] = _load_json(frozen_section44)
    config["formal_section44_protocol_source_path"] = str(frozen_section44.resolve())
    config["formal_section44_protocol_sha256"] = sha256(frozen_section44)
    config["top_level_entrypoint_id"] = entrypoint_id
    for dataset in config.get("datasets", []):
        suffix = str(dataset["id"]).split("stanford_vie_", 1)[-1]
        folder = preprocessing_root / suffix
        for attribute in dataset.get("attributes", []):
            property_id = str(attribute["property_id"])
            attribute["grid_path"] = str((folder / f"stanford_vie_{suffix}_{property_id}.grid").resolve())
            if property_id == "channel_occupancy":
                attribute["valid_mask_path"] = str((folder / f"stanford_vie_{suffix}_channel_label_valid.grid").resolve())
    runtime_path = run_dir / "configs" / "template_component_runtime.json"
    _write_json(runtime_path, config)
    return runtime_path


def _run_section44(python: Path, config: Path, mode: str, run_dir: Path) -> dict[str, Any]:
    component_mode = "formal" if mode == "formal" else "smoke"
    record = _run_command(
        [str(python), str(TEMPLATE_ROOT / "run_template_directional_experiment.py"), "--config", str(config), "--mode", component_mode, "--output-root", str(run_dir / "section_4_4")],
        run_dir / "logs" / "section44.log.txt",
        cwd=TEMPLATE_ROOT,
    )
    _require_command(record, "section44")
    gates = sorted((run_dir / "section_4_4").rglob("gate.json"), key=lambda path: path.stat().st_mtime)
    if not gates:
        raise RuntimeError("Section 4.4 did not produce gate.json")
    gate = _load_json(gates[-1])
    expected_status = "SMOKE_PASS" if component_mode == "smoke" else "DIAGNOSTIC"
    if gate.get("status") != expected_status or not all(bool(value) for value in gate.get("gates", {}).values()):
        raise RuntimeError(f"Section 4.4 gate did not pass: {gates[-1]}")
    return {
        "name": "section44",
        "passed": True,
        "details": str(gates[-1].resolve()),
        "gate_path": str(gates[-1].resolve()),
        "surrogate_count": int(gate["surrogate_count"]),
    }


def _run_direct_controls(python: Path, common_path: Path, common: Mapping[str, Any], mode: str, run_dir: Path) -> dict[str, Any]:
    source = _resolve(common_path, str(common["controlled_config"]))
    component_root = run_dir / "components" / "controlled_positive_controls"
    command = [str(python), str(CONTROLLED_ROOT / "run_direct_controlled_moment_experiment.py"), "--config", str(source), "--output-root", str(component_root)]
    if mode in {"smoke", "analysis-smoke"}:
        command.append("--quick")
    record = _run_command(command, run_dir / "logs" / "controlled_positive_controls.log.txt", cwd=CONTROLLED_ROOT)
    _require_command(record, "controlled_positive_controls")
    gates = sorted(component_root.rglob("gate.json"), key=lambda path: path.stat().st_mtime)
    if not gates:
        raise RuntimeError("Direct controlled component did not write gate.json")
    gate = _load_json(gates[-1])
    if not bool(gate.get("overall_pass")) or gate.get("status") != "VALIDATION_PASS":
        raise RuntimeError(f"Direct controlled component did not pass: {gates[-1]}")
    runtime = _load_json(gates[-1].parent / "run_config.json")
    return {
        "name": "controlled_positive_controls",
        "passed": True,
        "details": str(gates[-1].resolve()),
        "gate_path": str(gates[-1].resolve()),
        "random_seed_count": len(runtime["random_seeds"]),
        "windows_per_class": int(runtime["windows_per_class"]),
        "bootstrap_replicates": int(runtime["bootstrap_replicates"]),
    }


def _run_parity(executable: Path, grids: Mapping[str, GridData], spacing: tuple[float, float, float], run_dir: Path) -> dict[str, Any]:
    fixture = run_cpp_python_fixture_matrix_gate(executable, run_dir / "cpp_python_fixture_matrix", source_spacing=spacing)
    routes = run_uniform_three_route_audit(executable, run_dir / "uniform_three_route_audit", source_spacing=spacing)
    real = run_cpp_python_parity_gate(executable, grids["z45_56"], "z45_56", run_dir / "cpp_python_parity", source_spacing=spacing)
    if not fixture.passed or fixture.case_count != 60 or not routes.passed or routes.case_count != 6 or not real.passed or real.case_count != 12:
        raise RuntimeError("Direct protocol C++/Python parity gate did not complete its full fixed matrix")
    return {
        "name": "cpp_python_parity_gate",
        "passed": True,
        "details": f"fixture={fixture.report_path}; uniform_routes={routes.report_path}; real_input={real.report_path}",
    }


def _artifact_groups(artifacts: Iterable[CppFeatureMapArtifact]) -> dict[tuple[str, str, str, str, str, str], dict[int, CppFeatureMapArtifact]]:
    groups: dict[tuple[str, str, str, str, str, str], dict[int, CppFeatureMapArtifact]] = {}
    for artifact in artifacts:
        row = artifact.manifest_row
        key = (
            str(row["section"]), str(row["dataset_id"]), str(row["window"]), str(row["weighting"]),
            str(row["parameter_set"]), str(row["parameter_role"]),
        )
        if int(row["order"]) in groups.setdefault(key, {}):
            raise ValueError(f"Duplicate C++ direct map group/order: {key}/{row['order']}")
        groups[key][int(row["order"])] = artifact
    return groups


def _coordinate_values(grid: GridData, property_name: str) -> pd.DataFrame:
    frame = grid.coordinates.rename(columns={grid.property_name: property_name})
    return frame.loc[:, ["x", "y", "z", property_name]].copy()


def _window_high_fraction(
    grid: GridData,
    spec: CppFeatureMapSpec,
    coordinates: pd.DataFrame,
    *,
    threshold: float,
    variance_epsilon: float,
) -> pd.DataFrame:
    batch = extract_complete_windows(grid, tuple(spec.window_size), tuple(spec.move_step))
    configuration = WeightConfig.from_mapping(dict(spec.weight_config))
    weights = generate_weights(
        batch.values,
        configuration,
        window_size=tuple(spec.window_size),
        source_spacing=tuple(spec.source_spacing),
    )
    statistics = compute_weighted_statistics(
        batch, weights, threshold, variance_epsilon, cumulant_order=6
    ).statistics
    result = batch.coordinates.copy()
    result["high_fraction"] = statistics.high_fraction
    merged = coordinates.merge(result, on=["x", "y", "z"], how="left", validate="one_to_one")
    if merged["high_fraction"].isna().any():
        raise RuntimeError("Python high-fraction windows do not align to the C++ centres")
    return merged.loc[:, ["x", "y", "z", "high_fraction"]]


def _independent_raw_lambda4(
    grid: GridData,
    spec: CppFeatureMapSpec,
    coordinates: pd.DataFrame,
    *,
    variance_epsilon: float,
) -> pd.DataFrame:
    """Recalculate excess kurtosis from raw window values and raw weights.

    This deliberately does not derive a second expression from C++ ``kappa4``.
    It applies the configured weights to the original window samples, computes
    the weighted mean, second central moment, and fourth central moment in
    Python, then forms ``mu4 / mu2**2 - 3`` independently.
    """

    batch = extract_complete_windows(grid, tuple(spec.window_size), tuple(spec.move_step))
    configuration = WeightConfig.from_mapping(dict(spec.weight_config))
    weights = generate_weights(
        batch.values,
        configuration,
        window_size=tuple(spec.window_size),
        source_spacing=tuple(spec.source_spacing),
    )
    weight_sum = np.sum(weights, axis=1)
    safe_sum = np.where(weight_sum > 0.0, weight_sum, 1.0)
    mean = np.sum(weights * batch.values, axis=1) / safe_sum
    centered = batch.values - mean[:, None]
    variance = np.sum(weights * centered**2, axis=1) / safe_sum
    central_moment4 = np.sum(weights * centered**4, axis=1) / safe_sum
    valid = np.isfinite(variance) & (variance > float(variance_epsilon)) & (weight_sum > 0.0)
    excess = np.full(len(batch.coordinates), np.nan, dtype=float)
    excess[valid] = central_moment4[valid] / variance[valid] ** 2 - 3.0
    raw = batch.coordinates.copy()
    raw["python_raw_variance"] = variance
    raw["python_raw_central_moment4"] = central_moment4
    raw["python_raw_excess_kurtosis"] = excess
    raw["python_raw_variance_valid"] = valid
    merged = coordinates.merge(raw, on=["x", "y", "z"], how="left", validate="one_to_one")
    if merged["python_raw_variance_valid"].isna().any():
        raise RuntimeError("Independent raw lambda4 windows do not align to the C++ centres")
    return merged


def _independent_lambda4_audit(
    case: DirectCase,
    grid: GridData,
    spec: CppFeatureMapSpec,
    *,
    variance_epsilon: float,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> dict[str, Any]:
    raw = _independent_raw_lambda4(
        grid,
        spec,
        case.frame.loc[:, ["x", "y", "z", "lambda4", "analysis_variance_valid"]],
        variance_epsilon=variance_epsilon,
    )
    valid = (
        raw["analysis_variance_valid"].astype(bool).to_numpy()
        & raw["python_raw_variance_valid"].astype(bool).to_numpy()
        & np.isfinite(raw["lambda4"].to_numpy(dtype=float))
        & np.isfinite(raw["python_raw_excess_kurtosis"].to_numpy(dtype=float))
    )
    if not valid.any():
        raise RuntimeError(f"Independent raw lambda4 audit has no valid windows: {case.identity}")
    difference = (
        raw.loc[valid, "lambda4"].to_numpy(dtype=float)
        - raw.loc[valid, "python_raw_excess_kurtosis"].to_numpy(dtype=float)
    )
    reference = np.maximum(
        np.maximum(
            np.abs(raw.loc[valid, "lambda4"].to_numpy(dtype=float)),
            np.abs(raw.loc[valid, "python_raw_excess_kurtosis"].to_numpy(dtype=float)),
        ),
        1.0,
    )
    tolerance = float(absolute_tolerance) + float(relative_tolerance) * reference
    if absolute_tolerance < 0.0 or relative_tolerance < 0.0:
        raise ValueError("Independent lambda4 audit tolerances must be non-negative")
    return {
        **dict(zip(("section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role"), case.identity, strict=True)),
        "audit_method": "python_raw_weighted_central_moment4",
        "valid_count": int(valid.sum()),
        "max_abs_error": float(np.max(np.abs(difference))),
        "mean_abs_error": float(np.mean(np.abs(difference))),
        "max_relative_error": float(np.max(np.abs(difference) / reference)),
        "absolute_tolerance": float(absolute_tolerance),
        "relative_tolerance": float(relative_tolerance),
        "within_tolerance": bool(np.all(np.abs(difference) <= tolerance)),
        "finite": bool(np.isfinite(difference).all()),
    }


def _direct_case_from_group(
    group_key: tuple[str, str, str, str, str, str],
    maps: Mapping[int, CppFeatureMapArtifact],
    grid: GridData,
    references: Mapping[str, GridData],
    *,
    threshold: float,
    variance_epsilon: float,
    minimum_ess: float | None,
) -> tuple[DirectCase, dict[str, Any]]:
    if set(maps) != {1, 2, 3, 4, 5, 6}:
        raise ValueError(f"Primary direct case is incomplete: {group_key}")
    first = maps[1]
    frame = first.map_frame.loc[:, ["x", "y", "z"]].copy()
    for order in range(1, 7):
        artifact = maps[order]
        map_frame = artifact.map_frame.rename(columns={"value": f"kappa{order}"})
        frame = frame.merge(map_frame, on=["x", "y", "z"], how="inner", validate="one_to_one")
    sidecar = first.support_frame.loc[first.support_frame["map_included"].astype(bool), ["x", "y", "z", "valid", "effective_sample_size"]].copy()
    sidecar = sidecar.rename(columns={"valid": "cpp_valid"})
    frame = frame.merge(sidecar, on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["label"], "channel"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["label_valid"], "label_valid"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["occupancy"], "channel_occupancy"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(
        _window_high_fraction(grid, first.spec, frame.loc[:, ["x", "y", "z"]], threshold=threshold, variance_epsilon=variance_epsilon),
        on=["x", "y", "z"], how="inner", validate="one_to_one",
    )
    frame["channel"] = frame["channel"].astype(bool)
    frame["label_valid"] = frame["label_valid"].astype(bool)
    frame["cpp_valid"] = frame["cpp_valid"].astype(bool)
    standardized = derive_standardized_cumulants(frame, variance_epsilon)
    supported, diagnostics = build_analysis_support(
        standardized,
        required_features=FEATURES,
        variance_epsilon=variance_epsilon,
        minimum_ess=minimum_ess,
    )
    section, dataset_id, window, weighting, parameter_set, role = group_key
    return DirectCase(section, dataset_id, window, weighting, parameter_set, role, supported), diagnostics


def _sensitivity_case_from_group(
    group_key: tuple[str, str, str, str, str, str],
    maps: Mapping[int, CppFeatureMapArtifact],
    grid: GridData,
    references: Mapping[str, GridData],
    *,
    variance_epsilon: float,
    minimum_ess: float,
) -> DirectCase:
    if set(maps) != {2, 6}:
        raise ValueError(f"Sensitivity direct case is incomplete: {group_key}")
    first = maps[2]
    frame = first.map_frame.loc[:, ["x", "y", "z"]].copy()
    for order in (2, 6):
        frame = frame.merge(maps[order].map_frame.rename(columns={"value": f"kappa{order}"}), on=["x", "y", "z"], how="inner", validate="one_to_one")
    sidecar = first.support_frame.loc[first.support_frame["map_included"].astype(bool), ["x", "y", "z", "valid", "effective_sample_size"]].rename(columns={"valid": "cpp_valid"})
    frame = frame.merge(sidecar, on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["label"], "channel"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["label_valid"], "label_valid"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    frame = frame.merge(_coordinate_values(references["occupancy"], "channel_occupancy"), on=["x", "y", "z"], how="inner", validate="one_to_one")
    variance = frame["kappa2"].to_numpy(dtype=float)
    valid = np.isfinite(variance) & (variance > variance_epsilon)
    frame["variance"] = variance
    frame["lambda6"] = np.where(valid, frame["kappa6"].to_numpy(dtype=float) / variance**3, np.nan)
    frame["channel"] = frame["channel"].astype(bool)
    frame["label_valid"] = frame["label_valid"].astype(bool)
    frame["cpp_valid"] = frame["cpp_valid"].astype(bool)
    supported, _ = build_analysis_support(
        frame,
        required_features=("variance", "lambda6"),
        variance_epsilon=variance_epsilon,
        minimum_ess=minimum_ess,
    )
    section, dataset_id, window, weighting, parameter_set, role = group_key
    return DirectCase(section, dataset_id, window, weighting, parameter_set, role, supported)


def _analysis_support_scope(case: DirectCase) -> str:
    if case.section == "4.6":
        return (
            "section46_sensitivity_main_and_sidewings_common"
            if case.parameter_role == "sensitivity"
            else "section46_all_three_primary_common"
        )
    return "section45_individual"


def _write_support_masks(run_dir: Path, cases: Sequence[DirectCase], grids: Mapping[str, GridData]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for index, case in enumerate(cases):
        grid = grids[case.dataset_id]
        mask = np.zeros_like(grid.values, dtype=bool)
        x_index = {int(value): number for number, value in enumerate(grid.x)}
        y_index = {int(value): number for number, value in enumerate(grid.y)}
        z_index = {int(value): number for number, value in enumerate(grid.z)}
        selected = case.frame.loc[case.frame["analysis_valid"].astype(bool), ["x", "y", "z"]]
        for row in selected.itertuples(index=False):
            mask[z_index[int(row.z)], y_index[int(row.y)], x_index[int(row.x)]] = True
        name = f"{index:02d}_{case.section}_{case.dataset_id}_{case.window}_{case.weighting}_{case.parameter_role}".replace(".", "_")
        npz, metadata, details = write_mask_bundle(
            run_dir / "support_masks" / name,
            mask,
            coordinate_basis="run-local coarse z/y/x grid; C++ output centres",
            extra={"case_identity": list(case.identity)},
        )
        ess = case.frame.loc[case.frame["analysis_valid"].astype(bool), "effective_sample_size"].to_numpy(dtype=float)
        rows.append(
            {
                "row_type": "support_mask",
                "support_scope": _analysis_support_scope(case),
                "section": case.section,
                "dataset_id": case.dataset_id,
                "window": case.window,
                "weighting": case.weighting,
                "parameter_set": case.parameter_set,
                "parameter_role": case.parameter_role,
                "candidate_count": int(len(case.frame)),
                "analysis_valid_count": int(case.frame["analysis_valid"].sum()),
                "analysis_valid_rate": float(case.frame["analysis_valid"].mean()),
                "ess_mean": float(np.mean(ess)),
                "ess_min": float(np.min(ess)),
                "ess_over_n_mean": float(np.mean(ess / (int(case.window.split("x")[0]) * int(case.window.split("x")[1]) * int(case.window.split("x")[2])))),
                "mask_npz": str(npz.relative_to(run_dir)),
                "mask_metadata": str(metadata.relative_to(run_dir)),
                "mask_sha256": str(details["npz_sha256"]),
            }
        )
    return pd.DataFrame(rows)


def _coordinate_gaussian_profile_rows(spec: CppFeatureMapSpec) -> list[dict[str, Any]]:
    """Record the discrete central x/y/z profile of the primary spatial weight."""

    configuration = WeightConfig.from_mapping(dict(spec.weight_config))
    if configuration.active_mode() != "coordinate_gaussian":
        raise ValueError("A coordinate profile requires coordinate_gaussian weights")
    wx, wy, wz = (int(value) for value in spec.window_size)
    weights = coordinate_gaussian_window_weights(
        (wx, wy, wz), tuple(spec.source_spacing), float(configuration.bandwidth), float(configuration.min_weight)
    ).reshape(wz, wy, wx)
    centres = {"x": wx // 2, "y": wy // 2, "z": wz // 2}
    profiles = {
        "x": weights[centres["z"], centres["y"], :],
        "y": weights[centres["z"], :, centres["x"]],
        "z": weights[:, centres["y"], centres["x"]],
    }
    dimensions = {"x": wx, "y": wy, "z": wz}
    rows: list[dict[str, Any]] = []
    for axis, values in profiles.items():
        size = dimensions[axis]
        half_width = max(1.0, 0.5 * (size - 1))
        for index, value in enumerate(values):
            rows.append(
                {
                    "diagnostic_type": "coordinate_gaussian_central_axis_profile",
                    "section": "4.6",
                    "dataset_id": "window_definition",
                    "window": "x".join(str(value) for value in spec.window_size),
                    "weighting": "coordinate_gaussian",
                    "parameter_set": str(spec.parameter_set),
                    "axis": axis,
                    "axis_index": int(index),
                    "normalized_offset": float((index - 0.5 * (size - 1)) / half_width),
                    "weight": float(value),
                    "tail_weight_fraction": np.nan,
                    "analysis_window_count": np.nan,
                    "sample_count": np.nan,
                }
            )
    return rows


def _soft_upper_tail_weight_fraction(
    case: DirectCase,
    grid: GridData,
    spec: CppFeatureMapSpec,
    *,
    threshold: float,
) -> dict[str, Any]:
    """Calculate the allocated soft-weight share for porosity at or above the threshold."""

    configuration = WeightConfig.from_mapping(dict(spec.weight_config))
    if configuration.active_mode() != "soft_high_porosity":
        raise ValueError("Upper-tail allocation requires soft_high_porosity weights")
    batch = extract_complete_windows(grid, tuple(spec.window_size), tuple(spec.move_step))
    weights = generate_weights(
        batch.values,
        configuration,
        window_size=tuple(spec.window_size),
        source_spacing=tuple(spec.source_spacing),
    )
    support = case.frame.loc[:, ["x", "y", "z", "analysis_valid"]]
    aligned = batch.coordinates.merge(support, on=["x", "y", "z"], how="left", validate="one_to_one")
    if aligned["analysis_valid"].isna().any():
        raise RuntimeError("Soft-weight diagnostic windows do not align to the direct analysis support")
    selected = aligned["analysis_valid"].to_numpy(dtype=bool)
    selected_weights = weights[selected]
    selected_values = batch.values[selected]
    total_weight = float(np.sum(selected_weights))
    if not np.isfinite(total_weight) or total_weight <= 0.0:
        raise RuntimeError("Soft-weight diagnostic has no positive common-support weight")
    upper_weight = float(np.sum(selected_weights[selected_values >= float(threshold)]))
    return {
        "diagnostic_type": "soft_high_porosity_upper_tail_weight_fraction",
        "section": "4.6",
        "dataset_id": case.dataset_id,
        "window": case.window,
        "weighting": case.weighting,
        "parameter_set": case.parameter_set,
        "axis": "",
        "axis_index": np.nan,
        "normalized_offset": np.nan,
        "weight": np.nan,
        "tail_weight_fraction": upper_weight / total_weight,
        "analysis_window_count": int(selected.sum()),
        "sample_count": int(selected_weights.size),
    }


def _write_weighting_mechanism_diagnostics(
    run_dir: Path,
    section46_cases: Sequence[DirectCase],
    specs_by_case: Mapping[tuple[str, str, str, str, str, str], CppFeatureMapSpec],
    grids: Mapping[str, GridData],
    *,
    threshold: float,
) -> pd.DataFrame:
    """Persist prespecified spatial and value-weight mechanism diagnostics."""

    coordinate_cases = [case for case in section46_cases if case.weighting == "coordinate_gaussian"]
    soft_cases = [case for case in section46_cases if case.weighting == "soft_high_porosity"]
    if not coordinate_cases or not soft_cases:
        raise RuntimeError("Section 4.6 mechanism diagnostics require coordinate and soft primary cases")
    rows = _coordinate_gaussian_profile_rows(specs_by_case[coordinate_cases[0].identity])
    rows.extend(
        _soft_upper_tail_weight_fraction(
            case,
            grids[case.dataset_id],
            specs_by_case[case.identity],
            threshold=threshold,
        )
        for case in soft_cases
    )
    diagnostics = pd.DataFrame(rows)
    diagnostics.to_csv(run_dir / "weighting_mechanism_diagnostics.csv", index=False)
    return diagnostics


def _with_bootstrap_group(case: DirectCase, group: str) -> DirectCase:
    return DirectCase(
        case.section,
        case.dataset_id,
        case.window,
        case.weighting,
        case.parameter_set,
        case.parameter_role,
        case.frame,
        group,
    )


def _build_sensitivity_comparisons(
    primary_cases: Sequence[DirectCase],
    sensitivity_cases: Sequence[DirectCase],
) -> tuple[list[DirectCase], list[dict[str, Any]]]:
    """Pair each side-wing weight setting with its same-scheme main setting."""

    paired: list[DirectCase] = []
    diagnostics: list[dict[str, Any]] = []
    primary_by_key = {
        (case.dataset_id, case.window, case.weighting): case
        for case in primary_cases
        if case.weighting in {"coordinate_gaussian", "soft_high_porosity"}
    }
    for weighting in sorted({case.weighting for case in sensitivity_cases}):
        wings = [case for case in sensitivity_cases if case.weighting == weighting]
        for dataset_id in sorted({case.dataset_id for case in wings}):
            selected_wings = [case for case in wings if case.dataset_id == dataset_id]
            window = selected_wings[0].window
            main = primary_by_key.get((dataset_id, window, weighting))
            if main is None:
                raise ValueError(f"No primary main setting for sensitivity comparison: {dataset_id}/{window}/{weighting}")
            common, detail = apply_common_support([main, *selected_wings])
            group = f"sensitivity:{dataset_id}:{window}:{weighting}"
            paired.extend(_with_bootstrap_group(case, group) for case in common)
            diagnostics.append(
                {
                    "row_type": "common_support",
                    "support_scope": "section46_sensitivity_main_and_sidewings_common",
                    "section": "4.6",
                    "dataset_id": dataset_id,
                    "window": window,
                    "weighting": weighting,
                    "parameter_set": "main_and_sidewings_common",
                    "parameter_role": "sensitivity",
                    **detail,
                }
            )
    return paired, diagnostics


def _scoped_bootstrap(
    cases: Sequence[DirectCase],
    *,
    scope: str,
    features: Sequence[str],
    block_shape_xyz: tuple[int, int, int],
    replicates: int,
    seed: int,
    include_advantages: bool = True,
    include_weighting_effects: bool = True,
    include_parameter_sensitivity: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    started = time.perf_counter()
    result = joint_spatial_block_bootstrap(
        cases,
        features=features,
        block_shape_xyz=block_shape_xyz,
        replicates=replicates,
        seed=seed,
        include_advantages=include_advantages,
        include_weighting_effects=include_weighting_effects,
        include_parameter_sensitivity=include_parameter_sensitivity,
    )
    elapsed = time.perf_counter() - started
    summary = result.summary.copy()
    summary["bootstrap_scope"] = scope
    diagnostics = {
        "scope": scope,
        "elapsed_seconds": float(elapsed),
        "case_count": int(len(cases)),
        "replicates": int(replicates),
        "block_shape_xyz": "x".join(str(value) for value in block_shape_xyz),
        **result.diagnostics,
    }
    return summary, diagnostics


def _write_time_projection(
    run_dir: Path,
    manifest: pd.DataFrame,
    *,
    mode: str,
    matrix_elapsed_seconds: float,
    bootstrap_diagnostics: Mapping[str, Mapping[str, Any]],
    common: Mapping[str, Any],
    component_seconds: Mapping[str, float],
    formal_bootstrap_case_counts: Mapping[str, int],
) -> dict[str, Any]:
    """Project the complete formal run from a sealed analysis-smoke run.

    The projection separates map execution, each bootstrap scope, and the
    fixed workflow components.  Bootstrap work is scaled by both the formal
    replication count and the number of direct cases in the corresponding
    scope, so a three-subvolume formal run is not priced as one smoke case.
    """

    representative = manifest.loc[
        (manifest["dataset_id"] == "z45_56")
        & (manifest["window"] == "5x5x3")
        & (manifest["weighting"] == "uniform")
        & manifest["order"].isin([2, 5, 6])
    ].copy()
    if mode == "analysis-smoke" and len(representative) != 3:
        raise RuntimeError("Analysis smoke must include exactly the prespecified 3-map timing fixture")
    per_map = float(representative["cli_elapsed_seconds"].sum() / len(representative)) if not representative.empty else float("nan")
    baseline_bootstrap = {
        name: {
            "elapsed_seconds": float(value["elapsed_seconds"]),
            "replicates": int(value["replicates"]),
            "case_count": int(value["case_count"]),
            "seconds_per_replicate": float(value["elapsed_seconds"]) / int(value["replicates"]),
        }
        for name, value in bootstrap_diagnostics.items()
    }
    formal = common["bootstrap"]
    formal_replicates = {
        "primary": int(formal["primary_replicates"]),
        "scale": int(formal["scale_replicates"]),
        "large_block": int(formal["large_block_replicates"]),
        "parameter_sensitivity": int(formal["scale_replicates"]),
    }
    projected_bootstrap_by_scope: dict[str, dict[str, float | int]] = {}
    for scope, observed in baseline_bootstrap.items():
        observed_cases = int(observed["case_count"])
        formal_cases = int(formal_bootstrap_case_counts[scope])
        if observed_cases <= 0 or formal_cases <= 0:
            raise RuntimeError(f"Bootstrap time projection has an invalid case count for {scope}")
        case_scale = formal_cases / observed_cases
        projected_bootstrap_by_scope[scope] = {
            "analysis_smoke_case_count": observed_cases,
            "formal_case_count": formal_cases,
            "case_scale_factor": case_scale,
            "formal_replicates": int(formal_replicates[scope]),
            "estimated_seconds": float(observed["seconds_per_replicate"]) * int(formal_replicates[scope]) * case_scale,
        }
    projected_bootstrap = float(sum(float(item["estimated_seconds"]) for item in projected_bootstrap_by_scope.values()))
    map_scale = 102.0 / float(len(manifest))
    section44_scale = float(component_seconds.get("section44_formal_scale", 1.0))
    controls_scale = float(component_seconds.get("controlled_positive_controls_formal_scale", 1.0))
    fixed_component_projection = {
        "release_build_and_tests": float(component_seconds.get("release_build_and_tests", 0.0)),
        "preprocessing": float(component_seconds.get("preprocessing", 0.0)),
        "section44": float(component_seconds.get("section44", 0.0)) * section44_scale,
        "controlled_positive_controls": float(component_seconds.get("controlled_positive_controls", 0.0)) * controls_scale,
        "cpp_python_parity": float(component_seconds.get("cpp_python_parity", 0.0)),
        "direct_analysis_overhead": float(component_seconds.get("direct_analysis_overhead", 0.0)) * map_scale,
        "report_generation": float(component_seconds.get("report_generation", 0.0)) * map_scale,
        "evidence_sealing": float(component_seconds.get("evidence_sealing_estimate", 0.0)) * map_scale,
    }
    projected_fixed = float(sum(fixed_component_projection.values()))
    estimated_map_seconds = float(per_map * 102) if np.isfinite(per_map) else float("nan")
    estimated_total = estimated_map_seconds + projected_bootstrap + projected_fixed
    budget_seconds = float(common.get("formal_time_budget_seconds", 6.0 * 3600.0))
    payload = {
        "schema_version": "direct-time-projection.v2",
        "mode": mode,
        "matrix_elapsed_seconds": float(matrix_elapsed_seconds),
        "three_map_fixture": {
            "identity": "z45_56;5x5x3;uniform;orders=2,5,6",
            "map_count": int(len(representative)),
            "elapsed_seconds": float(representative["cli_elapsed_seconds"].sum()) if not representative.empty else float("nan"),
            "seconds_per_map": per_map,
        },
        "bootstrap_baseline": baseline_bootstrap,
        "bootstrap_projection_by_scope": projected_bootstrap_by_scope,
        "analysis_smoke_component_seconds": {
            key: float(value) for key, value in component_seconds.items()
            if not key.endswith("_formal_scale")
        },
        "formal_component_scale_factors": {
            "map_or_map_derived_output": map_scale,
            "section44_surrogate_work": section44_scale,
            "controlled_positive_control_work": controls_scale,
        },
        "formal_fixed_component_projection_seconds": fixed_component_projection,
        "formal_projection": {
            "map_count": 102,
            "estimated_map_seconds": estimated_map_seconds,
            "estimated_bootstrap_seconds": projected_bootstrap,
            "estimated_fixed_workflow_seconds": projected_fixed,
            "estimated_total_seconds": estimated_total,
            "budget_seconds": budget_seconds,
            "time_budget_pass": bool(np.isfinite(estimated_total) and estimated_total <= budget_seconds),
        },
    }
    _write_json(run_dir / "time_projection.json", payload)
    return payload


def _run_direct_analysis(
    specs: Sequence[CppFeatureMapSpec],
    *,
    executable: Path,
    grids: Mapping[str, GridData],
    references: Mapping[str, Mapping[str, GridData]],
    common: Mapping[str, Any],
    run_dir: Path,
    mode: str,
) -> dict[str, Any]:
    if mode == "formal":
        counts = validate_direct_cpp_feature_map_protocol(specs)
    else:
        counts = {
            "map_count": len(specs),
            "primary_count": int(sum(spec.parameter_role == "primary" for spec in specs)),
            "sensitivity_count": int(sum(spec.parameter_role == "sensitivity" for spec in specs)),
        }
    matrix_started = time.perf_counter()
    result = run_cpp_feature_matrix(executable, specs, run_dir, timeout_seconds=900.0, manifest_name="cpp_feature_map_manifest.csv")
    matrix_elapsed_seconds = time.perf_counter() - matrix_started
    if mode == "smoke":
        return {"map_counts": counts, "cpp_manifest": str(result.manifest_path), "smoke_map_only": True}
    grouped = _artifact_groups(result.artifacts)
    primary_cases: list[DirectCase] = []
    sensitivity_cases: list[DirectCase] = []
    support_rows: list[dict[str, Any]] = []
    equivalence_rows: list[dict[str, Any]] = []
    primary_specs_by_case: dict[tuple[str, str, str, str, str, str], CppFeatureMapSpec] = {}
    for group_key, maps in sorted(grouped.items()):
        section, dataset_id, _, weighting, _, role = group_key
        if role == "primary":
            minimum_ess = None if weighting == "uniform" else float(common["minimum_ess_order_6"])
            case, diagnostics = _direct_case_from_group(
                group_key, maps, grids[dataset_id], references[dataset_id],
                threshold=float(common["high_porosity_threshold"]), variance_epsilon=float(common["variance_epsilon"]), minimum_ess=minimum_ess,
            )
            primary_cases.append(case)
            primary_specs_by_case[case.identity] = maps[1].spec
            support_rows.append(
                {
                    "row_type": "raw_support",
                    "support_scope": "individual_weighting_before_common_support",
                    **dict(zip(("section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role"), case.identity, strict=True)),
                    **diagnostics,
                }
            )
            equivalence_rows.append(
                _independent_lambda4_audit(
                    case,
                    grids[dataset_id],
                    maps[1].spec,
                    variance_epsilon=float(common["variance_epsilon"]),
                    absolute_tolerance=float(common["lambda4_independent_audit_absolute_tolerance"]),
                    relative_tolerance=float(common["lambda4_independent_audit_relative_tolerance"]),
                )
            )
        else:
            sensitivity_cases.append(
                _sensitivity_case_from_group(
                    group_key, maps, grids[dataset_id], references[dataset_id],
                    variance_epsilon=float(common["variance_epsilon"]), minimum_ess=float(common["minimum_ess_order_6"]),
                )
            )
    section45_cases = [case for case in primary_cases if case.section == "4.5"]
    section46_weighted = [case for case in primary_cases if case.section == "4.6"]
    uniform_by_dataset = {(case.dataset_id, case.window): case for case in section45_cases if case.weighting == "uniform"}
    section46_cases: list[DirectCase] = []
    for dataset_id in sorted({case.dataset_id for case in section46_weighted}):
        uniform = uniform_by_dataset[(dataset_id, "5x5x3")]
        companions = [uniform, *[case for case in section46_weighted if case.dataset_id == dataset_id]]
        common_cases, common_details = apply_common_support(companions)
        uniform_common = common_cases[0]
        section46_cases.append(DirectCase("4.6", uniform_common.dataset_id, uniform_common.window, uniform_common.weighting, uniform_common.parameter_set, uniform_common.parameter_role, uniform_common.frame))
        section46_cases.extend(common_cases[1:])
        support_rows.append(
            {
                "row_type": "common_support",
                "support_scope": "section46_all_three_primary_common",
                "section": "4.6",
                "dataset_id": dataset_id,
                "window": "5x5x3",
                "weighting": "all_three",
                "parameter_set": "common",
                "parameter_role": "primary",
                **common_details,
            }
        )
    all_primary = [*section45_cases, *section46_cases]
    effects, feature_summary = compute_direct_effects(all_primary, FEATURES)
    effects = annotate_signed_direction_consistency(effects)
    advantages = compute_single_map_advantages(effects)
    weighting_effects = compute_weighting_effects(effects, advantages)
    sensitivity_pairs, sensitivity_support = _build_sensitivity_comparisons(primary_cases, sensitivity_cases)
    support_rows.extend(sensitivity_support)
    sensitivity_effects, _ = compute_direct_effects(sensitivity_pairs, ("variance", "lambda6"))
    sensitivity_effects = annotate_signed_direction_consistency(sensitivity_effects, expected_dataset_count=1)
    parameter_sensitivity = compute_parameter_sensitivity_effects(sensitivity_effects)

    bootstrap_config = common["bootstrap"]
    replicates = (
        {
            "primary": ANALYSIS_SMOKE_BOOTSTRAP_REPLICATES,
            "scale": ANALYSIS_SMOKE_BOOTSTRAP_REPLICATES,
            "large_block": ANALYSIS_SMOKE_BOOTSTRAP_REPLICATES,
            "parameter_sensitivity": ANALYSIS_SMOKE_BOOTSTRAP_REPLICATES,
        }
        if mode == "analysis-smoke"
        else {
            "primary": int(bootstrap_config["primary_replicates"]),
            "scale": int(bootstrap_config["scale_replicates"]),
            "large_block": int(bootstrap_config["large_block_replicates"]),
            "parameter_sensitivity": int(bootstrap_config["scale_replicates"]),
        }
    )
    seed = int(bootstrap_config["seed"])
    primary_window_cases = [case for case in section45_cases if case.window == "5x5x3"]
    primary_window_cases.extend(section46_cases)
    scale_cases = [case for case in section45_cases if case.window in {"3x3x3", "7x7x3"}]
    primary_summary, primary_diagnostics = _scoped_bootstrap(
        primary_window_cases,
        scope="primary",
        features=FEATURES,
        block_shape_xyz=tuple(int(value) for value in common["block_shape_primary_xyz"]),
        replicates=replicates["primary"],
        seed=seed,
    )
    scale_summary, scale_diagnostics = _scoped_bootstrap(
        scale_cases,
        scope="scale",
        features=FEATURES,
        block_shape_xyz=tuple(int(value) for value in common["block_shape_primary_xyz"]),
        replicates=replicates["scale"],
        seed=seed + 1,
        include_weighting_effects=False,
    )
    large_summary, large_diagnostics = _scoped_bootstrap(
        primary_window_cases,
        scope="large_block",
        features=FEATURES,
        block_shape_xyz=tuple(int(value) for value in common["block_shape_large_xyz"]),
        replicates=replicates["large_block"],
        seed=seed + 2,
    )
    sensitivity_summary, sensitivity_diagnostics = _scoped_bootstrap(
        sensitivity_pairs,
        scope="parameter_sensitivity",
        features=("variance", "lambda6"),
        block_shape_xyz=tuple(int(value) for value in common["block_shape_primary_xyz"]),
        replicates=replicates["parameter_sensitivity"],
        seed=seed + 3,
        include_advantages=False,
        include_weighting_effects=False,
        include_parameter_sensitivity=True,
    )
    bootstrap_summary = pd.concat(
        (primary_summary, scale_summary, large_summary, sensitivity_summary),
        ignore_index=True,
        sort=False,
    )
    pd.DataFrame(equivalence_rows).to_csv(run_dir / "lambda4_independent_audit.csv", index=False)
    if not all(row["finite"] and row["within_tolerance"] for row in equivalence_rows):
        raise RuntimeError("Independent raw-window lambda4/excess-kurtosis audit failed")
    effects.to_csv(run_dir / "direct_effects.csv", index=False)
    feature_summary.to_csv(run_dir / "direct_feature_summary.csv", index=False)
    advantages.to_csv(run_dir / "single_map_advantages.csv", index=False)
    weighting_effects.to_csv(run_dir / "weighting_effects.csv", index=False)
    parameter_sensitivity.to_csv(run_dir / "parameter_sensitivity_contrasts.csv", index=False)
    bootstrap_summary.to_csv(run_dir / "bootstrap_summary.csv", index=False)
    # The side-wing pairing includes its primary reference cases.  The primary
    # maps are already represented in all_primary, so only side-wing masks are
    # added here.  This preserves one ESS/N observation per formal map family.
    sensitivity_masks = [case for case in sensitivity_pairs if case.parameter_role == "sensitivity"]
    mask_summary = _write_support_masks(run_dir, [*all_primary, *sensitivity_masks], grids)
    support = pd.concat((pd.DataFrame(support_rows), mask_summary), ignore_index=True, sort=False)
    support.to_csv(run_dir / "support_summary.csv", index=False)
    _write_weighting_mechanism_diagnostics(
        run_dir,
        section46_cases,
        primary_specs_by_case,
        grids,
        threshold=float(common["high_porosity_threshold"]),
    )
    # Section 4.6 is persisted after common-support intersection.  The report
    # therefore draws the identical coordinate support used for W/G estimates.
    values = pd.concat(
        [
            case.frame.assign(
                section=case.section,
                dataset_id=case.dataset_id,
                window=case.window,
                weighting=case.weighting,
                parameter_set=case.parameter_set,
                parameter_role=case.parameter_role,
                analysis_support_scope=_analysis_support_scope(case),
            )
            for case in [*section45_cases, *section46_cases]
        ],
        ignore_index=True,
    )
    values.to_csv(run_dir / "direct_analysis_values.csv", index=False)
    bootstrap_diagnostics = {
        "primary": primary_diagnostics,
        "scale": scale_diagnostics,
        "large_block": large_diagnostics,
        "parameter_sensitivity": sensitivity_diagnostics,
    }
    return {
        "map_counts": counts,
        "cpp_manifest": str(result.manifest_path),
        "matrix_elapsed_seconds": float(matrix_elapsed_seconds),
        "lambda4_equivalence_max_abs_error": float(max(row["max_abs_error"] for row in equivalence_rows)),
        "bootstrap_diagnostics": bootstrap_diagnostics,
        "primary_case_count": len(all_primary),
        "sensitivity_case_count": len(sensitivity_pairs),
    }


def _write_report(run_dir: Path, *, mode: str, status: str, direct: Mapping[str, Any], checks: Sequence[Mapping[str, Any]]) -> None:
    lines = [
        "# Direct-Statistics Sections 4.4-4.6 Run Report",
        "",
        f"Mode: `{mode}`",
        f"Status: `{status}`",
        "",
        "## Direct Protocol",
        "",
        "- Sections 4.5 and 4.6 use direct map statistics only; no classifier, cross-validation, learned threshold, fold, or OOF artifact is generated.",
        f"- C++ feature maps: `{direct['map_counts']['map_count']}`",
        f"- Primary maps: `{direct['map_counts']['primary_count']}`",
        f"- Sensitivity maps: `{direct['map_counts']['sensitivity_count']}`",
    ]
    if "lambda4_equivalence_max_abs_error" in direct:
        lines.append(f"- lambda4/excess-kurtosis maximum absolute error: `{direct['lambda4_equivalence_max_abs_error']:.3e}`")
    lines.extend(["", "## Checks", ""])
    lines.extend(f"- `{item['name']}`: {'PASS' if item['passed'] else 'FAIL'}" for item in checks)
    lines.extend(["", "## Core Artifacts", ""])
    lines.extend(f"- `{name}`" for name in DIRECT_REQUIRED_ARTIFACTS if name not in {"gate.json", "evidence_manifest.json", "sha256_manifest.csv"})
    (run_dir / "RUN_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_results_index(run_dir: Path, *, mode: str, checks: Sequence[Mapping[str, Any]]) -> None:
    lines = [
        "# Direct-Statistics Result Index",
        "",
        f"Mode: `{mode}`",
        "",
        "## Outputs",
        "",
        "- `preprocessing/preprocessing_manifest.json`",
        "- `section_4_4/`",
        "- `components/controlled_positive_controls/`",
        "- `cpp_feature_map_manifest.csv`",
        "- `direct_feature_summary.csv`, `direct_effects.csv`, `single_map_advantages.csv`, `weighting_effects.csv`",
        "- `bootstrap_summary.csv`, `support_summary.csv`, `common_support_manifest.csv`",
        "- `lambda4_independent_audit.csv`, `parameter_sensitivity_contrasts.csv`, `weighting_mechanism_diagnostics.csv`",
        "- `report_artifact_contract.csv`",
        "- `RUN_REPORT.md`, `gate.json`, `evidence_manifest.json`, `sha256_manifest.csv`",
        "",
        "## Components",
        "",
    ]
    lines.extend(f"- `{item['name']}`: {item['details']}" for item in checks)
    (run_dir / "RESULTS_INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _forbidden_output_paths(run_dir: Path) -> list[str]:
    """Return only explicit historical model-protocol artifacts."""

    forbidden: list[str] = []
    for path in run_dir.rglob("*"):
        if not path.is_file():
            continue
        name = path.name.lower()
        relative = path.relative_to(run_dir)
        # Section 4.4 retains its established surrogate omnibus/ Holm table.
        # Its file name predates the direct protocol but contains no model,
        # fold, OOF, coefficient, or learned-threshold output.
        if relative.parts[:1] == ("section_4_4",) and name == "multiplicity_decisions.csv":
            continue
        if name in FORBIDDEN_OUTPUT_ARTIFACTS or any(name.startswith(prefix) and name.endswith(".csv") for prefix in FORBIDDEN_OUTPUT_PREFIXES):
            forbidden.append(relative.as_posix())
    return sorted(forbidden)


def _direct_output_contract(run_dir: Path, *, mode: str, sealed: bool) -> dict[str, Any]:
    required = list(DIRECT_REQUIRED_ARTIFACTS)
    if mode == "smoke":
        required = ["cpp_feature_map_manifest.csv", "RUN_REPORT.md", "RESULTS_INDEX.md", "gate.json", "evidence_manifest.json", "sha256_manifest.csv"]
    if not sealed:
        required = [name for name in required if name not in {"gate.json", "sha256_manifest.csv"}]
    missing = [name for name in required if not (run_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"Direct output contract is missing artifacts: {missing}")
    forbidden = _forbidden_output_paths(run_dir)
    if forbidden:
        raise RuntimeError(f"Direct output contract contains prohibited model artifacts: {forbidden}")
    manifest = pd.read_csv(run_dir / "cpp_feature_map_manifest.csv")
    centering_columns = {"order", "algorithm", "cumulant_centering", "cumulant_centering_version"}
    missing_centering = centering_columns.difference(manifest.columns)
    if missing_centering:
        raise RuntimeError(
            "Direct map manifest is missing cumulant-centering provenance: "
            f"{sorted(missing_centering)}"
        )
    if set(manifest["cumulant_centering_version"].astype(str)) != {CUMULANT_CENTERING_VERSION}:
        raise RuntimeError("Direct map manifest contains an unsupported cumulant-centering version")
    expected_centering = np.where(
        manifest["order"].astype(int) == 1,
        "none",
        np.where(
            manifest["algorithm"].astype(str) == "Weighted moment",
            "weighted_window_mean",
            "window_mean",
        ),
    )
    if not np.array_equal(manifest["cumulant_centering"].astype(str).to_numpy(), expected_centering):
        raise RuntimeError("Direct map manifest cumulant-centering modes do not match order and algorithm")
    if mode == "formal":
        if len(manifest) != 102:
            raise RuntimeError(f"Direct formal map manifest must contain 102 maps, found {len(manifest)}")
        full_key = ["section", "dataset_id", "window", "weighting_id", "parameter_set", "parameter_role", "order", "algorithm"]
        if manifest.duplicated(full_key).any():
            raise RuntimeError("Direct formal map manifest contains duplicate full identities")
    elif mode == "analysis-smoke":
        if len(manifest) != 34:
            raise RuntimeError(f"Direct analysis smoke must contain 34 maps, found {len(manifest)}")
    return {
        "required_artifact_count": len(required),
        "forbidden_artifact_count": len(forbidden),
        "map_count": int(len(manifest)),
        "sealed": bool(sealed),
    }


def _common_support_manifest(run_dir: Path) -> None:
    support = pd.read_csv(run_dir / "support_summary.csv")
    required = {"row_type", "support_scope", "section", "parameter_role"}
    missing = required.difference(support.columns)
    if missing:
        raise RuntimeError(f"Support summary is missing explicit support metadata: {sorted(missing)}")
    selected = support.loc[
        (support["section"].astype(str) == "4.6")
        & support["support_scope"].isin(
            ["section46_all_three_primary_common", "section46_sensitivity_main_and_sidewings_common"]
        )
        & support["row_type"].isin(["common_support", "support_mask"])
    ].copy()
    if selected.empty or not (selected["row_type"] == "common_support").any():
        raise RuntimeError("Common-support manifest has no common-support summary rows")
    if not (selected["row_type"] == "support_mask").any():
        raise RuntimeError("Common-support manifest has no coordinate-level support masks")
    selected.to_csv(run_dir / "common_support_manifest.csv", index=False)


def _source_files(config_sources: Mapping[str, Path]) -> list[Path]:
    files = [
        Path(__file__), R15_ROOT / "cli_runner.py", R15_ROOT / "direct_statistical_evaluation.py", R15_ROOT / "report_direct_results.py",
        R15_ROOT / "weighting.py", R15_ROOT / "statistics.py", R15_ROOT / "cpp_python_parity_gate.py",
        EXPERIMENTS_ROOT / "verify_direct_evidence.py", EXPERIMENTS_ROOT / "configure_paths.py",
        PREPROCESS_ROOT / "prepare_stanford_vie.py", TEMPLATE_ROOT / "formal_runner.py", TEMPLATE_ROOT / "run_template_directional_experiment.py",
        CONTROLLED_ROOT / "direct_evaluation.py", CONTROLLED_ROOT / "run_direct_controlled_moment_experiment.py",
    ]
    files.extend(config_sources.values())
    return files


def _analysis_smoke_specs(specs: Sequence[CppFeatureMapSpec]) -> tuple[CppFeatureMapSpec, ...]:
    """Return the reduced matrix that still exercises every analysis branch."""

    local = [spec for spec in specs if spec.dataset_id == "z45_56"]
    primary_specs = [spec for spec in local if spec.parameter_role == "primary"]
    sensitivity_by_mode: dict[str, list[CppFeatureMapSpec]] = {}
    for spec in local:
        if spec.parameter_role == "sensitivity":
            sensitivity_by_mode.setdefault(str(spec.weighting_id), []).append(spec)
    selected_sensitivity: list[CppFeatureMapSpec] = []
    for mode, candidates in sorted(sensitivity_by_mode.items()):
        parameter_set = min(str(spec.parameter_set) for spec in candidates)
        selected_sensitivity.extend(spec for spec in candidates if str(spec.parameter_set) == parameter_set)
    selected = tuple([*primary_specs, *selected_sensitivity])
    primary = sum(spec.parameter_role == "primary" for spec in selected)
    sensitivity = sum(spec.parameter_role == "sensitivity" for spec in selected)
    if (len(selected), primary, sensitivity) != (34, 30, 4):
        raise RuntimeError(
            "Direct analysis smoke must contain 30 primary maps and four side-wing maps"
        )
    return selected


def _formal_bootstrap_case_counts(specs: Sequence[CppFeatureMapSpec]) -> dict[str, int]:
    """Count direct-map case families that participate in each formal scope."""

    groups = {
        (
            spec.section,
            spec.dataset_id,
            "x".join(str(value) for value in spec.window_size),
            spec.weighting_id,
            spec.parameter_set,
            spec.parameter_role,
        )
        for spec in specs
    }
    section45 = {group for group in groups if group[0] == "4.5" and group[5] == "primary"}
    section46 = {group for group in groups if group[0] == "4.6" and group[5] == "primary"}
    if not section45 or not section46:
        raise RuntimeError("Formal bootstrap case count requires Section 4.5 and Section 4.6 primary groups")
    primary_window = "5x5x3"
    section45_primary = {group for group in section45 if group[2] == primary_window}
    section45_scale = {group for group in section45 if group[2] != primary_window}
    datasets = {group[1] for group in section46}
    # The 4.6 uniform map is inherited from the Section 4.5 primary window and
    # re-evaluated on the all-three common-support mask.
    section46_with_uniform = len(section46) + len(datasets)
    sensitivity = {group for group in groups if group[5] == "sensitivity"}
    sensitivity_main = {(dataset, window, weighting) for _, dataset, window, weighting, _, _ in sensitivity}
    counts = {
        "primary": len(section45_primary) + section46_with_uniform,
        "scale": len(section45_scale),
        "large_block": len(section45_primary) + section46_with_uniform,
        "parameter_sensitivity": len(sensitivity) + len(sensitivity_main),
    }
    if any(value <= 0 for value in counts.values()):
        raise RuntimeError(f"Formal bootstrap case counts are incomplete: {counts}")
    return counts


def _estimate_evidence_sealing_seconds(run_dir: Path) -> dict[str, float | int]:
    """Measure one archive hash pass and project the four sealing hash passes."""

    started = time.perf_counter()
    paths = [path for path in run_dir.rglob("*") if path.is_file() and path.name != "sha256_manifest.csv"]
    byte_count = 0
    for path in paths:
        byte_count += int(path.stat().st_size)
        sha256(path)
    one_pass = time.perf_counter() - started
    return {
        "hash_inventory_file_count": len(paths),
        "hash_inventory_bytes": byte_count,
        "hash_inventory_seconds": float(one_pass),
        # seal_evidence_bundle builds and verifies the manifest twice.
        "sealing_estimate_seconds": float(max(1.0, 4.0 * one_pass)),
    }


def run(mode: str, output_root: Path, config_root: Path = CONFIG_ROOT) -> Path:
    if mode not in {"smoke", "analysis-smoke", "formal"}:
        raise ValueError("mode must be smoke, analysis-smoke, or formal")
    config_root = config_root.resolve()
    common_source = config_root / "direct_common_protocol.json"
    common = _load_json(common_source)
    config_sources = {
        "direct_common_protocol": common_source,
        "direct_section45": config_root / "direct_section45.json",
        "direct_section46": config_root / "direct_section46.json",
        "formal_section44": config_root / "formal_section44.json",
        "preprocessing_config": _resolve(common_source, str(common["preprocessing_config"])),
        "template_component_config": _resolve(common_source, str(common["template_component_config"])),
        "controlled_component_config": _resolve(common_source, str(common["controlled_config"])),
    }
    if any(not path.is_file() for path in config_sources.values()):
        raise FileNotFoundError("A direct-statistics protocol configuration is missing")
    common = _load_json(common_source)
    section45 = _load_json(config_sources["direct_section45"])
    section46 = _load_json(config_sources["direct_section46"])
    section44_protocol = _load_json(config_sources["formal_section44"])
    controlled_protocol = _load_json(config_sources["controlled_component_config"])
    run_dir = _create_run(output_root)
    frozen = _copy_frozen_configs(run_dir, config_sources)
    configuration_hash = _combined_config_hash(config_sources)
    entrypoint_id = f"direct-44-46-{uuid.uuid4()}-{configuration_hash[:12]}"
    run_manifest = {
        "entrypoint_id": entrypoint_id,
        "entrypoint_path": str(Path(__file__).resolve()),
        "mode": mode,
        "frozen_config_sha256": configuration_hash,
        "frozen_config_files": {name: {"path": str(path.relative_to(run_dir)), "sha256": sha256(path)} for name, path in frozen.items()},
    }
    _write_json(run_dir / "run_manifest.json", run_manifest)
    python = Path(sys.executable).resolve()
    component_seconds: dict[str, float] = {}
    started = time.perf_counter()
    checks = _build_release_and_tests(python, common, run_dir, common_path=common_source)
    component_seconds["release_build_and_tests"] = time.perf_counter() - started
    started = time.perf_counter()
    preprocessing = _run_preprocessing(python, common_source, common, run_dir)
    component_seconds["preprocessing"] = time.perf_counter() - started
    checks.append({key: value for key, value in preprocessing.items() if key != "manifest"})
    grids, references, _ = _load_direct_grids(common_source, common, run_dir / "preprocessing")
    template_config = _runtime_template_config(common_source, common, frozen["formal_section44"], run_dir / "preprocessing", run_dir, entrypoint_id)
    started = time.perf_counter()
    section44_check = _run_section44(python, template_config, mode, run_dir)
    component_seconds["section44"] = time.perf_counter() - started
    checks.append(section44_check)
    started = time.perf_counter()
    controls_check = _run_direct_controls(python, common_source, common, mode, run_dir)
    component_seconds["controlled_positive_controls"] = time.perf_counter() - started
    checks.append(controls_check)
    executable = _resolve(common_source, str(common["movewindow_cli"]))
    spacing = tuple(float(value) for value in common["source_spacing_xyz_m"])
    started = time.perf_counter()
    checks.append(_run_parity(executable, grids, spacing, run_dir))
    component_seconds["cpp_python_parity"] = time.perf_counter() - started
    formal_specs = build_direct_cpp_feature_map_specs(
        grids,
        section45_windows=[tuple(int(value) for value in window) for window in section45["windows"]],
        section46_primary_window=tuple(int(value) for value in section46["primary_window"]),
        uniform_weight=dict(section45["weightings"][0]),
        nonuniform_primary_weights=[dict(weight) for weight in section46["primary_weightings"] if str(weight["mode"]) != "uniform"],
        sensitivity_weights=[dict(weight) for weight in section46["sensitivity_weightings"]],
        orders=tuple(int(value) for value in section45["orders"]),
        sensitivity_orders=tuple(int(value) for value in section46["sensitivity_orders"]),
        sensitivity_dataset_id=str(section46["sensitivity_dataset_id"]),
        move_step=tuple(int(value) for value in common["move_step"]),
        source_spacing=spacing,
    )
    formal_bootstrap_case_counts = _formal_bootstrap_case_counts(formal_specs)
    specs = formal_specs
    if mode == "smoke":
        specs = tuple(spec for spec in specs if spec.dataset_id == "z45_56" and tuple(spec.window_size) == (5, 5, 3) and spec.weight_config.get("mode") == "uniform" and int(spec.order) in {2, 5, 6})
    elif mode == "analysis-smoke":
        # One fixed subvolume retains all three Section 4.5 scales, both
        # nonuniform primary schemes, and one prespecified side wing for each
        # weighting family.  It is the smallest matrix that exercises common
        # support, paired side-wing contrasts, bootstrap aggregation, tables,
        # and figures end to end.
        specs = _analysis_smoke_specs(specs)
    started = time.perf_counter()
    direct = _run_direct_analysis(specs, executable=executable, grids=grids, references=references, common=common, run_dir=run_dir, mode=mode)
    direct_elapsed = time.perf_counter() - started
    bootstrap_elapsed = sum(float(item["elapsed_seconds"]) for item in direct.get("bootstrap_diagnostics", {}).values())
    component_seconds["direct_analysis_overhead"] = max(
        0.0, direct_elapsed - float(direct.get("matrix_elapsed_seconds", 0.0)) - bootstrap_elapsed
    )
    checks.append({"name": "direct_cpp_map_matrix", "passed": True, "details": str(direct["cpp_manifest"])})
    if mode in {"formal", "analysis-smoke"}:
        _common_support_manifest(run_dir)
        started = time.perf_counter()
        report = generate_direct_report(run_dir)
        component_seconds["report_generation"] = time.perf_counter() - started
        checks.append({"name": "direct_report_artifacts", "passed": True, "details": str(report["artifact_contract"])})
        status = "PASS" if mode == "formal" else "ANALYSIS_SMOKE_PASS"
    else:
        (run_dir / "common_support_manifest.csv").write_text("section,dataset_id,window,weighting\n", encoding="utf-8")
        status = "SMOKE_PASS"
        component_seconds["report_generation"] = 0.0
    observed_surrogates = int(section44_check["surrogate_count"])
    if observed_surrogates <= 0:
        raise RuntimeError("Section 4.4 timing projection requires a positive observed surrogate count")
    component_seconds["section44_formal_scale"] = float(section44_protocol["formal_surrogates_per_dataset_property"]) / observed_surrogates
    control_seed_scale = len(controlled_protocol["random_seeds"]) / int(controls_check["random_seed_count"])
    control_window_scale = int(controlled_protocol["train_windows_per_class"]) * 2 / int(controls_check["windows_per_class"])
    control_bootstrap_scale = (
        len(controlled_protocol["random_seeds"]) / int(controls_check["random_seed_count"])
        * int(controlled_protocol["bootstrap_replicates"]) / int(controls_check["bootstrap_replicates"])
    )
    component_seconds["controlled_positive_controls_formal_scale"] = max(
        control_seed_scale * control_window_scale,
        control_bootstrap_scale,
    )
    sealing_estimate = _estimate_evidence_sealing_seconds(run_dir)
    component_seconds["evidence_sealing_estimate"] = float(sealing_estimate["sealing_estimate_seconds"])
    if mode in {"formal", "analysis-smoke"}:
        direct["time_projection"] = _write_time_projection(
            run_dir,
            pd.read_csv(run_dir / "cpp_feature_map_manifest.csv"),
            mode=mode,
            matrix_elapsed_seconds=float(direct["matrix_elapsed_seconds"]),
            bootstrap_diagnostics=direct["bootstrap_diagnostics"],
            common=common,
            component_seconds=component_seconds,
            formal_bootstrap_case_counts=formal_bootstrap_case_counts,
        )
        direct["evidence_sealing_estimate"] = sealing_estimate
        budget_pass = bool(direct["time_projection"]["formal_projection"]["time_budget_pass"])
        checks.append(
            {
                "name": "formal_time_budget_projection",
                "passed": budget_pass,
                "details": str(run_dir / "time_projection.json"),
            }
        )
        if not budget_pass:
            raise RuntimeError("Projected complete formal workflow exceeds the configured time budget")
    _write_report(run_dir, mode=mode, status=status, direct=direct, checks=checks)
    _write_results_index(run_dir, mode=mode, checks=checks)
    environment = collect_environment(executable, _source_files(config_sources), int(common["bootstrap"]["seed"]))
    evidence_manifest = {
        "schema_version": "direct-sections-44-46.v1",
        "mode": mode,
        "entrypoint_id": entrypoint_id,
        "run_manifest": {"path": "run_manifest.json", "sha256": sha256(run_dir / "run_manifest.json")},
        "preprocessing_manifest": {"path": "preprocessing/preprocessing_manifest.json", "sha256": sha256(run_dir / "preprocessing" / "preprocessing_manifest.json")},
        "environment": environment,
        "checks": checks,
        "direct": direct,
    }
    _write_json(run_dir / "evidence_manifest.json", evidence_manifest)
    # The pre-seal contract covers all produced analysis artifacts.  The final
    # sealed pass below additionally verifies the gate and hash manifest.
    contract = _direct_output_contract(run_dir, mode=mode, sealed=False)
    checks.append({"name": "direct_output_contract", "passed": True, "details": json.dumps(contract, sort_keys=True)})
    evidence_manifest["checks"] = checks
    evidence_manifest["output_contract"] = contract
    _write_json(run_dir / "evidence_manifest.json", evidence_manifest)
    seal_evidence_bundle(
        run_dir,
        checks,
        status,
        metadata={"mode": mode, "entrypoint_id": entrypoint_id, "evaluation_mode": "direct_statistics", "map_counts": direct["map_counts"]},
    )
    _direct_output_contract(run_dir, mode=mode, sealed=True)
    return run_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "analysis-smoke", "formal"), required=True)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "output" / "direct_4_4_4_6_20261011")
    parser.add_argument("--config-root", type=Path, default=CONFIG_ROOT, help="Configuration directory generated by configure_paths.py.")
    arguments = parser.parse_args(argv)
    result = run(str(arguments.mode), arguments.output_root, arguments.config_root)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
