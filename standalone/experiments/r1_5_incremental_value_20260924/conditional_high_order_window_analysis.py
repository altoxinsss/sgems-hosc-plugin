"""Conditional empirical-window analysis for the Section 4.5 uniform maps.

This is an intentionally small, post-processing-only analysis.  It consumes a
completed direct-statistics run, recovers the original 5x5x3 porosity windows,
and asks what local marginal-distribution differences remain after lower-order
statistics are conditionally accounted for.  It does not call either C++ CLI,
does not use the facies labels, and does not fit a predictor.

The script never writes into ``run_dir``.  Keeping derivative evidence in a
separate directory preserves the manifest of the completed formal run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from direct_statistical_evaluation import build_spatial_blocks
from grid_io import read_grid
from statistics import local_statistics_from_values
from window_support import extract_complete_windows


SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENTS_ROOT = SCRIPT_PATH.parents[1]
COORDINATES = ("x", "y", "z")
WINDOW = (5, 5, 3)
WINDOW_ID = "5x5x3"
UNIFORM_WINDOW_SAMPLE_SIZE = int(np.prod(WINDOW))
ILLUSTRATIVE_PAIR_NEIGHBORS = 8
LOWER_CONTROLS_5 = ("mean", "variance", "high_fraction", "lambda3", "lambda4")
LOWER_CONTROLS_6 = (*LOWER_CONTROLS_5, "lambda5")
TAIL_THRESHOLDS = (-1.5, 1.5)
ECDF_GRID = np.linspace(-9.0, 9.0, 361)
CONTROL_DISPLAY_NAMES = {"variance": "log_variance"}

# These three associations are fixed before inspecting the conditional results.
# They follow the roles of the standardized cumulants: lambda5 describes an
# asymmetric fifth-order component, while lambda6 describes sixth-order shape
# after lambda5 has been accounted for.  The complete tail-metric matrix is
# still exported for audit rather than hiding the remaining descriptors.
PRIMARY_ASSOCIATION_METRICS: Mapping[str, tuple[str, ...]] = {
    "lambda5": ("tail_imbalance",),
    "lambda6": ("outer_tail_probability", "outer_span_q05_q95"),
}


@dataclass(frozen=True)
class PreparedDataset:
    """One coordinate-aligned uniform-window family."""

    dataset_id: str
    frame: pd.DataFrame
    values: np.ndarray
    standardized_values: np.ndarray


@dataclass(frozen=True)
class TargetResult:
    """Rank-projection and residual-group result for one cumulant order."""

    dataset_id: str
    target: str
    controls: tuple[str, ...]
    rank: np.ndarray
    residual: np.ndarray
    residual_q20: float
    residual_q80: float
    group: np.ndarray
    r_squared: float
    coefficients: np.ndarray


@dataclass(frozen=True)
class RankOLSResult:
    """Rank transform and rank-OLS residual for one response variable."""

    rank: np.ndarray
    residual: np.ndarray
    r_squared: float
    coefficients: np.ndarray


@dataclass(frozen=True)
class AssociationBlockAggregates:
    """Blockwise sufficient statistics for a fixed-rank conditional association.

    ``z`` contains the target rank followed by the ranks of all empirical
    distribution descriptors.  Recombining these matrices in a block bootstrap
    refits both nuisance projections without re-reading the source windows.
    """

    dataset_id: str
    target: str
    controls: tuple[str, ...]
    metric_names: tuple[str, ...]
    xtx: np.ndarray
    xtz: np.ndarray
    ztz: np.ndarray


@dataclass(frozen=True)
class FrozenRunInputs:
    """Immutable formal-run inputs used by this derivative analysis."""

    common: dict[str, Any]
    common_config_path: Path
    common_config_sha256: str
    run_manifest_path: Path
    run_manifest_sha256: str
    gate_path: Path
    gate_sha256: str
    evidence_manifest_path: Path
    evidence_manifest_sha256: str
    archived_hash_manifest_path: Path
    archived_hash_manifest_sha256: str
    preprocessing_manifest_path: Path
    preprocessing_manifest_sha256: str
    porosity_paths: Mapping[str, Path]
    porosity_hashes: Mapping[str, str]
    reference_variance: Mapping[str, float]
    variance_epsilon: Mapping[str, float]


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_dump(path: Path, value: Mapping[str, Any] | list[Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _archived_hash_index(run_directory: Path) -> tuple[Path, dict[str, str]]:
    """Load the final formal-run SHA-256 inventory without trusting its paths."""

    manifest_path = run_directory / "sha256_manifest.csv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Formal run lacks sha256_manifest.csv: {manifest_path}")
    frame = pd.read_csv(manifest_path, usecols=["path", "sha256"])
    if frame.empty or frame["path"].duplicated().any():
        raise ValueError("Formal SHA-256 inventory is empty or contains duplicate paths")
    values = {
        str(path).replace("\\", "/"): str(digest).lower()
        for path, digest in zip(frame["path"], frame["sha256"], strict=True)
    }
    if not all(len(digest) == 64 for digest in values.values()):
        raise ValueError("Formal SHA-256 inventory contains an invalid digest")
    return manifest_path, values


def _verify_archived_file(path: Path, run_directory: Path, inventory: Mapping[str, str], description: str) -> str:
    """Verify one formal-run artifact against its final inventory."""

    resolved = path.resolve()
    _require_within(resolved, run_directory, description)
    if not resolved.is_file():
        raise FileNotFoundError(f"{description} is absent: {resolved}")
    key = str(resolved.relative_to(run_directory)).replace("\\", "/")
    expected = inventory.get(key)
    if expected is None:
        raise ValueError(f"{description} is absent from the formal SHA-256 inventory: {key}")
    observed = sha256(resolved)
    if observed != expected:
        raise RuntimeError(f"{description} hash does not match the formal SHA-256 inventory")
    return observed


def _verify_formal_run(run_directory: Path) -> dict[str, Any]:
    """Require a completed, archived direct-statistics run before post-processing."""

    hash_manifest_path, inventory = _archived_hash_index(run_directory)
    run_manifest_path = run_directory / "run_manifest.json"
    gate_path = run_directory / "gate.json"
    evidence_manifest_path = run_directory / "evidence_manifest.json"
    values_path = run_directory / "direct_analysis_values.csv"
    run_manifest_hash = _verify_archived_file(run_manifest_path, run_directory, inventory, "run_manifest.json")
    gate_hash = _verify_archived_file(gate_path, run_directory, inventory, "gate.json")
    evidence_hash = _verify_archived_file(evidence_manifest_path, run_directory, inventory, "evidence_manifest.json")
    values_hash = _verify_archived_file(values_path, run_directory, inventory, "direct_analysis_values.csv")
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    evidence = json.loads(evidence_manifest_path.read_text(encoding="utf-8"))
    if str(run_manifest.get("mode")) != "formal":
        raise ValueError("Conditional analysis requires a formal source run")
    if str(gate.get("mode")) != "formal" or str(gate.get("status")) != "PASS":
        raise ValueError("Conditional analysis requires a PASS formal gate")
    checks = gate.get("checks", [])
    if not checks or not all(isinstance(item, dict) and bool(item.get("passed")) for item in checks):
        raise ValueError("Formal gate has missing or failed checks")
    if str(evidence.get("schema_version")) != "direct-sections-44-46.v1":
        raise ValueError("Formal evidence manifest has an incompatible schema")
    if str(evidence.get("mode")) != "formal":
        raise ValueError("Formal evidence manifest does not describe a formal run")
    for name, expected_hash in (("run_manifest", run_manifest_hash), ("preprocessing_manifest", None)):
        entry = evidence.get(name, {})
        if not isinstance(entry, dict):
            raise ValueError(f"Formal evidence manifest lacks {name} metadata")
        declared = str(entry.get("sha256", "")).lower()
        if name == "run_manifest" and declared != expected_hash:
            raise RuntimeError("Evidence manifest run_manifest hash mismatch")
    try:
        section45_entry = run_manifest["frozen_config_files"]["direct_section45"]
        section45_relative = Path(str(section45_entry["path"]))
        section45_expected_hash = str(section45_entry["sha256"]).lower()
    except (KeyError, TypeError) as error:
        raise ValueError("run_manifest.json lacks frozen direct_section45 metadata") from error
    if section45_relative.is_absolute():
        raise ValueError("Frozen Section 4.5 configuration path must be run-relative")
    section45_path = (run_directory / section45_relative).resolve()
    _require_within(section45_path, run_directory, "Frozen Section 4.5 configuration")
    if _verify_archived_file(section45_path, run_directory, inventory, "Frozen Section 4.5 configuration") != section45_expected_hash:
        raise RuntimeError("Frozen Section 4.5 configuration hash mismatch")
    section45 = json.loads(section45_path.read_text(encoding="utf-8"))
    if (
        tuple(int(value) for value in section45.get("primary_window", ())) != WINDOW
        or [str(item.get("id")) for item in section45.get("weightings", [])] != ["uniform"]
        or tuple(int(value) for value in section45.get("orders", ())) != (1, 2, 3, 4, 5, 6)
    ):
        raise ValueError("Frozen Section 4.5 configuration is incompatible with the conditional analysis")
    return {
        "hash_manifest_path": hash_manifest_path,
        "hash_manifest_sha256": sha256(hash_manifest_path),
        "inventory": inventory,
        "run_manifest_path": run_manifest_path,
        "run_manifest_sha256": run_manifest_hash,
        "gate_path": gate_path,
        "gate_sha256": gate_hash,
        "evidence_manifest_path": evidence_manifest_path,
        "evidence_manifest_sha256": evidence_hash,
        "direct_analysis_values_path": values_path,
        "direct_analysis_values_sha256": values_hash,
        "direct_section45_path": section45_path,
        "direct_section45_sha256": section45_expected_hash,
        "evidence": evidence,
    }


def _bool_values(series: pd.Series) -> np.ndarray:
    """Handle native booleans and CSV spellings without truthy-string errors."""

    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(dtype=bool)
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "yes"}).to_numpy(dtype=bool)


def _require_within(path: Path, parent: Path, description: str) -> None:
    """Reject paths outside a frozen formal-run subtree."""

    try:
        path.relative_to(parent)
    except ValueError as error:
        raise ValueError(f"{description} is outside the frozen formal run: {path}") from error


def _load_frozen_common_protocol(
    run_directory: Path,
    supplied_path: Path | None,
) -> tuple[dict[str, Any], Path, Path, str]:
    """Load the protocol copied into the formal run and verify its frozen hash.

    A user-supplied configuration is accepted only when its content hash is
    identical to the copy named by ``run_manifest.json``.  Its data paths are
    intentionally ignored: all grids are recovered from the preprocessing
    snapshot under the supplied formal run.
    """

    run_manifest_path = run_directory / "run_manifest.json"
    if not run_manifest_path.is_file():
        raise FileNotFoundError(f"Formal run lacks run_manifest.json: {run_manifest_path}")
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    try:
        frozen = run_manifest["frozen_config_files"]["direct_common_protocol"]
        relative_path = Path(str(frozen["path"]))
        expected_hash = str(frozen["sha256"]).lower()
    except (KeyError, TypeError) as error:
        raise ValueError("run_manifest.json lacks frozen direct_common_protocol metadata") from error
    if relative_path.is_absolute():
        raise ValueError("Frozen common-protocol path must be relative to the formal run")
    frozen_path = (run_directory / relative_path).resolve()
    _require_within(frozen_path, run_directory, "Frozen common-protocol path")
    if not frozen_path.is_file():
        raise FileNotFoundError(f"Frozen common protocol is absent: {frozen_path}")
    frozen_hash = sha256(frozen_path)
    if frozen_hash != expected_hash:
        raise RuntimeError("Frozen common-protocol hash does not match run_manifest.json")
    selected_path = frozen_path if supplied_path is None else supplied_path.resolve()
    if not selected_path.is_file():
        raise FileNotFoundError(f"Supplied common protocol is absent: {selected_path}")
    if sha256(selected_path) != expected_hash:
        raise ValueError("Supplied common protocol differs from the formal run's frozen protocol")
    payload = json.loads(selected_path.read_text(encoding="utf-8"))
    if tuple(int(value) for value in payload.get("move_step", ())) != (1, 1, 1):
        raise ValueError("Conditional analysis requires the formal unit-step moving-window support")
    high_threshold = float(payload.get("high_porosity_threshold", np.nan))
    if not np.isfinite(high_threshold):
        raise ValueError("Frozen common protocol has a non-finite high-porosity threshold")
    dataset_ids: set[str] = set()
    for item in payload.get("datasets", []):
        identifier = str(item.get("id", "")).strip()
        if not identifier or identifier in dataset_ids:
            raise ValueError("Common protocol contains absent or duplicate dataset ids")
        dataset_ids.add(identifier)
    if not dataset_ids:
        raise ValueError("Common protocol contains no dataset ids")
    return payload, frozen_path, run_manifest_path, expected_hash


def _load_run_local_porosity_snapshot(
    run_directory: Path,
    expected_dataset_ids: Iterable[str],
    archived_inventory: Mapping[str, str],
) -> tuple[Path, dict[str, Path], dict[str, str]]:
    """Resolve and verify run-local porosity grids from preprocessing evidence."""

    preprocessing_root = (run_directory / "preprocessing").resolve()
    manifest_path = preprocessing_root / "preprocessing_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Formal run lacks preprocessing manifest: {manifest_path}")
    _verify_archived_file(manifest_path, run_directory, archived_inventory, "Preprocessing manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    schema = str(manifest.get("schema_version", ""))
    if not schema.startswith("stanford-vie-preprocessing."):
        raise ValueError(f"Incompatible preprocessing manifest schema: {schema!r}")
    paths: dict[str, Path] = {}
    hashes: dict[str, str] = {}
    for subvolume in manifest.get("subvolumes", []):
        identifier = str(subvolume.get("subvolume_id", "")).strip()
        outputs = [
            item for item in subvolume.get("outputs", [])
            if str(item.get("property_name", "")).strip().lower() == "porosity"
        ]
        if not identifier:
            raise ValueError("Preprocessing manifest contains an unnamed subvolume")
        if len(outputs) != 1:
            raise ValueError(f"Preprocessing manifest requires exactly one porosity grid for {identifier}")
        output = outputs[0]
        path = Path(str(output.get("path", ""))).resolve()
        _require_within(path, preprocessing_root, f"Porosity grid for {identifier}")
        if not path.is_file():
            raise FileNotFoundError(f"Run-local porosity grid is absent: {path}")
        expected_hash = str(output.get("sha256", "")).lower()
        observed_hash = _verify_archived_file(path, run_directory, archived_inventory, f"Porosity grid for {identifier}")
        if len(expected_hash) != 64 or observed_hash != expected_hash:
            raise RuntimeError(f"Run-local porosity-grid hash verification failed for {identifier}")
        if identifier in paths:
            raise ValueError(f"Preprocessing manifest contains duplicate subvolume id: {identifier}")
        paths[identifier] = path
        hashes[identifier] = observed_hash
    expected = set(expected_dataset_ids)
    if set(paths) != expected:
        raise ValueError(
            "Frozen preprocessing datasets do not match the frozen common protocol: "
            f"expected={sorted(expected)}, observed={sorted(paths)}"
        )
    return manifest_path, paths, hashes


def _derive_variance_thresholds(porosity_paths: Mapping[str, Path]) -> tuple[dict[str, float], dict[str, float]]:
    """Reproduce the formal per-dataset lambda eligibility threshold exactly."""

    reference_variance: dict[str, float] = {}
    variance_epsilon: dict[str, float] = {}
    for identifier, path in porosity_paths.items():
        grid = read_grid(path)
        values = np.asarray(grid.values[grid.present], dtype=float)
        if values.size == 0 or not np.isfinite(values).all():
            raise ValueError(f"Run-local porosity grid has no finite values: {identifier}")
        reference = float(np.var(values, ddof=0))
        if not np.isfinite(reference) or reference < 0.0:
            raise ValueError(f"Run-local porosity grid has invalid reference variance: {identifier}")
        reference_variance[identifier] = reference
        variance_epsilon[identifier] = max(1.0e-12, 1.0e-8 * reference)
    return reference_variance, variance_epsilon


def _load_frozen_inputs(
    run_directory: Path,
    supplied_path: Path | None,
    formal_bundle: Mapping[str, Any],
) -> FrozenRunInputs:
    common, common_path, run_manifest_path, common_hash = _load_frozen_common_protocol(run_directory, supplied_path)
    if _verify_archived_file(common_path, run_directory, formal_bundle["inventory"], "Frozen common protocol") != common_hash:
        raise RuntimeError("Frozen common protocol is inconsistent with the formal SHA-256 inventory")
    dataset_ids = [str(item["id"]) for item in common["datasets"]]
    preprocessing_manifest_path, porosity_paths, porosity_hashes = _load_run_local_porosity_snapshot(
        run_directory, dataset_ids, formal_bundle["inventory"]
    )
    evidence_preprocessing = formal_bundle["evidence"].get("preprocessing_manifest", {})
    if (
        not isinstance(evidence_preprocessing, dict)
        or str(evidence_preprocessing.get("sha256", "")).lower() != sha256(preprocessing_manifest_path)
    ):
        raise RuntimeError("Evidence manifest preprocessing hash mismatch")
    reference_variance, variance_epsilon = _derive_variance_thresholds(porosity_paths)
    return FrozenRunInputs(
        common=common,
        common_config_path=common_path,
        common_config_sha256=common_hash,
        run_manifest_path=run_manifest_path,
        run_manifest_sha256=str(formal_bundle["run_manifest_sha256"]),
        gate_path=Path(formal_bundle["gate_path"]),
        gate_sha256=str(formal_bundle["gate_sha256"]),
        evidence_manifest_path=Path(formal_bundle["evidence_manifest_path"]),
        evidence_manifest_sha256=str(formal_bundle["evidence_manifest_sha256"]),
        archived_hash_manifest_path=Path(formal_bundle["hash_manifest_path"]),
        archived_hash_manifest_sha256=str(formal_bundle["hash_manifest_sha256"]),
        preprocessing_manifest_path=preprocessing_manifest_path,
        preprocessing_manifest_sha256=sha256(preprocessing_manifest_path),
        porosity_paths=porosity_paths,
        porosity_hashes=porosity_hashes,
        reference_variance=reference_variance,
        variance_epsilon=variance_epsilon,
    )


def _formal_block_shape(common: Mapping[str, Any], supplied_shape: Iterable[int] | None) -> tuple[int, int, int]:
    """Use the formal primary block geometry, rejecting incompatible overrides."""

    frozen = tuple(int(value) for value in common.get("block_shape_primary_xyz", ()))
    if len(frozen) != 3 or any(value <= 0 for value in frozen):
        raise ValueError("Frozen common protocol has an invalid primary block shape")
    if supplied_shape is None:
        return frozen
    provided = tuple(int(value) for value in supplied_shape)
    if provided != frozen:
        raise ValueError(
            f"Conditional analysis block shape must match frozen formal protocol: expected {frozen}, got {provided}"
        )
    return frozen


def _select_numeric_support(
    csv_path: Path,
    variance_epsilon: Mapping[str, float],
    reference_variance: Mapping[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Section 4.5 numerical support without using any facies fields.

    ``analysis_valid`` is deliberately not used: it includes ``label_valid`` in
    the formal direct-statistics workflow.  This postprocessor needs only the
    coordinate-aligned numerical support of the uniform cumulant maps.
    """

    required_columns = [
        *COORDINATES,
        "kappa1", "kappa2", "kappa3", "kappa4", "kappa5", "kappa6",
        "cpp_valid", "effective_sample_size", "mean", "variance", "high_fraction",
        "lambda3", "lambda4", "lambda5", "lambda6", "analysis_variance_valid",
        "section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role",
    ]
    accepted_datasets = set(variance_epsilon)
    if set(reference_variance) != accepted_datasets:
        raise ValueError("Reference-variance and epsilon dataset identities differ")
    chunks: list[pd.DataFrame] = []
    diagnostics: dict[str, dict[str, Any]] = {
        dataset_id: {
            "dataset_id": dataset_id,
            "structural_rows": 0,
            "nonfinite_excluded": 0,
            "cpp_invalid_excluded": 0,
            "effective_sample_size_not_75_excluded": 0,
            "variance_at_or_below_epsilon_excluded": 0,
            "formal_analysis_variance_flag_disagreement": 0,
        }
        for dataset_id in sorted(accepted_datasets)
    }
    for chunk in pd.read_csv(csv_path, usecols=required_columns, chunksize=100_000):
        section = pd.to_numeric(chunk["section"], errors="coerce")
        mask = (
            section.eq(4.5)
            & chunk["dataset_id"].astype(str).isin(accepted_datasets)
            & chunk["window"].astype(str).eq(WINDOW_ID)
            & chunk["weighting"].astype(str).eq("uniform")
            & chunk["parameter_set"].astype(str).eq("uniform")
            & chunk["parameter_role"].astype(str).eq("primary")
        )
        if not mask.any():
            continue
        selected = chunk.loc[mask].copy()
        selected["dataset_id"] = selected["dataset_id"].astype(str)
        numeric = [
            "kappa1", "kappa2", "kappa3", "kappa4", "kappa5", "kappa6",
            "effective_sample_size", "mean", "variance", "high_fraction",
            "lambda3", "lambda4", "lambda5", "lambda6",
        ]
        for name in numeric:
            selected[name] = pd.to_numeric(selected[name], errors="coerce")
        finite = np.isfinite(selected[numeric].to_numpy(dtype=float)).all(axis=1)
        cpp_valid = _bool_values(selected["cpp_valid"])
        exact_ess = np.isclose(
            selected["effective_sample_size"].to_numpy(dtype=float),
            float(UNIFORM_WINDOW_SAMPLE_SIZE),
            rtol=0.0,
            atol=1.0e-12,
        )
        epsilon = selected["dataset_id"].map(variance_epsilon).to_numpy(dtype=float)
        variance_valid = selected["kappa2"].to_numpy(dtype=float) > epsilon
        formal_variance_valid = _bool_values(selected["analysis_variance_valid"])
        for dataset_id, row_indices in selected.groupby("dataset_id", sort=False).groups.items():
            index = np.asarray(list(row_indices), dtype=int)
            # GroupBy preserves the chunk's original labels, so align through a
            # positional mask instead of treating those labels as array offsets.
            group = selected.loc[index]
            group_finite = np.isfinite(group[numeric].to_numpy(dtype=float)).all(axis=1)
            group_cpp = _bool_values(group["cpp_valid"])
            group_ess = np.isclose(
                group["effective_sample_size"].to_numpy(dtype=float),
                float(UNIFORM_WINDOW_SAMPLE_SIZE), rtol=0.0, atol=1.0e-12,
            )
            group_epsilon = float(variance_epsilon[str(dataset_id)])
            group_variance = group["kappa2"].to_numpy(dtype=float) > group_epsilon
            group_formal = _bool_values(group["analysis_variance_valid"])
            row = diagnostics[str(dataset_id)]
            row["structural_rows"] += int(len(group))
            row["nonfinite_excluded"] += int((~group_finite).sum())
            row["cpp_invalid_excluded"] += int((group_finite & ~group_cpp).sum())
            row["effective_sample_size_not_75_excluded"] += int((group_finite & group_cpp & ~group_ess).sum())
            row["variance_at_or_below_epsilon_excluded"] += int((group_finite & group_cpp & group_ess & ~group_variance).sum())
            row["formal_analysis_variance_flag_disagreement"] += int(
                (group_finite & group_cpp & group_ess & (group_formal != group_variance)).sum()
            )
        valid = finite & cpp_valid & exact_ess & variance_valid
        retained = selected.loc[valid].copy()
        if not retained.empty:
            chunks.append(retained)
    if not chunks:
        raise ValueError("No Section 4.5 uniform 5x5x3 numerical support was found")
    result = pd.concat(chunks, ignore_index=True)
    if set(result["dataset_id"].astype(str)) != accepted_datasets:
        missing = accepted_datasets.difference(result["dataset_id"].astype(str))
        raise ValueError(f"Direct values omit requested datasets: {sorted(missing)}")
    if result.duplicated(["dataset_id", *COORDINATES]).any():
        raise ValueError("Selected direct values contain duplicate dataset/coordinate records")
    result = result.sort_values(["dataset_id", "z", "y", "x"], kind="stable").reset_index(drop=True)
    support_rows: list[dict[str, Any]] = []
    for dataset_id in sorted(accepted_datasets):
        retained = result.loc[result["dataset_id"].astype(str) == dataset_id]
        if retained.empty:
            raise ValueError(f"No independently eligible numerical windows remain for {dataset_id}")
        row = diagnostics[dataset_id]
        row.update(
            {
                "reference_variance": float(reference_variance[dataset_id]),
                "variance_epsilon": float(variance_epsilon[dataset_id]),
                "required_effective_sample_size": UNIFORM_WINDOW_SAMPLE_SIZE,
                "numerical_windows": int(len(retained)),
                "window_samples": UNIFORM_WINDOW_SAMPLE_SIZE,
                "minimum_ess": float(retained["effective_sample_size"].min()),
                "maximum_ess": float(retained["effective_sample_size"].max()),
                "minimum_kappa2_retained": float(retained["kappa2"].min()),
                "minimum_variance_retained": float(retained["variance"].min()),
                "maximum_variance_retained": float(retained["variance"].max()),
                "label_fields_used": False,
            }
        )
        support_rows.append(row)
    return result, pd.DataFrame(support_rows)


def _recover_dataset_windows(frame: pd.DataFrame, porosity_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Recover the exact 75 raw values behind every selected map coordinate."""

    grid = read_grid(porosity_path)
    batch = extract_complete_windows(grid, WINDOW, (1, 1, 1))
    complete_coordinates = batch.coordinates
    lookup = pd.MultiIndex.from_frame(complete_coordinates.loc[:, list(COORDINATES)])
    requested = pd.MultiIndex.from_frame(frame.loc[:, list(COORDINATES)])
    positions = lookup.get_indexer(requested)
    if np.any(positions < 0):
        missing = frame.loc[positions < 0, list(COORDINATES)].head().to_dict("records")
        raise ValueError(f"Direct map centres are absent from recovered complete windows: {missing}")
    return batch.values[positions], positions


def _relative_error(observed: np.ndarray, expected: np.ndarray) -> np.ndarray:
    scale = np.maximum(np.maximum(np.abs(observed), np.abs(expected)), 1.0e-15)
    return np.abs(observed - expected) / scale


def parity_audit(
    prepared: PreparedDataset,
    *,
    high_threshold: float,
    variance_epsilon: float,
    sample_size: int,
    rtol: float,
    atol: float,
) -> pd.DataFrame:
    """Recompute deterministic windows and fail on an estimator/path mismatch."""

    count = len(prepared.frame)
    indices = np.unique(np.linspace(0, count - 1, min(count, sample_size), dtype=int))
    recomputed = local_statistics_from_values(
        prepared.values[indices],
        high_threshold=high_threshold,
        variance_epsilon=variance_epsilon,
        max_order=6,
    ).as_dict()
    direct = prepared.frame.iloc[indices]
    source_columns = {
        "mean": "mean",
        "variance": "variance",
        "lambda3": "skewness",
        "lambda4": "excess_kurtosis",
        "lambda5": "standardized_kappa5",
        "lambda6": "standardized_kappa6",
    }
    rows: list[dict[str, Any]] = []
    for direct_name, computed_name in source_columns.items():
        actual = direct[direct_name].to_numpy(dtype=float)
        expected = np.asarray(recomputed[computed_name], dtype=float)
        difference = np.abs(actual - expected)
        relative = _relative_error(actual, expected)
        passed = bool(np.allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=False))
        rows.append(
            {
                "dataset_id": prepared.dataset_id,
                "field": direct_name,
                "sampled_windows": int(len(indices)),
                "max_abs_error": float(difference.max()),
                "max_relative_error": float(relative.max()),
                "rtol": float(rtol),
                "atol": float(atol),
                "passed": passed,
            }
        )
    report = pd.DataFrame(rows)
    if not report["passed"].all():
        failures = report.loc[~report["passed"], ["dataset_id", "field", "max_abs_error", "max_relative_error"]]
        raise RuntimeError(f"Recovered-window parity audit failed: {failures.to_dict('records')}")
    return report


def _standardize_windows(values: np.ndarray, variance_epsilon: float) -> np.ndarray:
    mean = np.mean(values, axis=1)
    variance = np.mean((values - mean[:, None]) ** 2, axis=1)
    if np.any(~np.isfinite(variance)) or np.any(variance <= float(variance_epsilon)):
        raise ValueError("Recovered numerical support contains a window below the formal variance threshold")
    return (values - mean[:, None]) / np.sqrt(variance)[:, None]


def _tail_metrics(standardized_values: np.ndarray) -> pd.DataFrame:
    lower, upper = TAIL_THRESHOLDS
    quantiles = np.quantile(standardized_values, [0.05, 0.10, 0.90, 0.95], axis=1).T
    p_lower = np.mean(standardized_values <= lower, axis=1)
    p_upper = np.mean(standardized_values >= upper, axis=1)
    return pd.DataFrame(
        {
            "lower_tail_probability": p_lower,
            "upper_tail_probability": p_upper,
            "outer_tail_probability": p_lower + p_upper,
            "tail_imbalance": p_upper - p_lower,
            "q05": quantiles[:, 0],
            "q10": quantiles[:, 1],
            "q90": quantiles[:, 2],
            "q95": quantiles[:, 3],
            "outer_span_q05_q95": quantiles[:, 3] - quantiles[:, 0],
        }
    )


def _unit_rank(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Rank projection requires finite values")
    return rankdata(values, method="average") / float(len(values) + 1)


def _control_label(name: str) -> str:
    """Return the public label for a rank-adjustment variable."""

    return CONTROL_DISPLAY_NAMES.get(name, name)


def _rank_ols_projection(response: np.ndarray, design: np.ndarray) -> RankOLSResult:
    """Project a fixed full-subvolume rank response on fixed rank controls."""

    response_rank = _unit_rank(response)
    coefficient, _, _, _ = np.linalg.lstsq(design, response_rank, rcond=None)
    residual = response_rank - design @ coefficient
    total = float(np.sum((response_rank - response_rank.mean()) ** 2))
    r_squared = 1.0 - float(np.sum(residual**2)) / total if total > 0.0 else float("nan")
    return RankOLSResult(
        rank=response_rank,
        residual=residual,
        r_squared=float(r_squared),
        coefficients=coefficient,
    )


def _rank_design(rank_controls: pd.DataFrame, controls: tuple[str, ...]) -> np.ndarray:
    """Build an intercept-plus-control rank design with an explicit order."""

    missing = [name for name in controls if f"rank_{name}" not in rank_controls]
    if missing:
        raise ValueError(f"Rank controls are missing columns: {missing}")
    return np.column_stack(
        (
            np.ones(len(rank_controls), dtype=float),
            *(rank_controls[f"rank_{name}"].to_numpy(dtype=float) for name in controls),
        )
    )


def _rank_response_projection(
    response: np.ndarray,
    rank_controls: pd.DataFrame,
    controls: tuple[str, ...],
) -> RankOLSResult:
    """Rank-residualize one response using the specified lower-order controls."""

    return _rank_ols_projection(np.asarray(response, dtype=float), _rank_design(rank_controls, controls))


def _rank_projection(
    frame: pd.DataFrame,
    target: str,
    controls: tuple[str, ...],
    *,
    rank_controls: pd.DataFrame | None = None,
) -> TargetResult:
    """Fit a deterministic rank-OLS projection for a cumulant target."""

    controls_frame = _rank_control_frame(frame) if rank_controls is None else rank_controls
    fit = _rank_response_projection(frame[target].to_numpy(dtype=float), controls_frame, controls)
    residual = fit.residual
    q20, q80 = np.quantile(residual, (0.20, 0.80))
    group = np.full(len(frame), "middle", dtype=object)
    group[residual <= q20] = "low_q20"
    group[residual >= q80] = "high_q80"
    return TargetResult(
        dataset_id=str(frame["dataset_id"].iloc[0]),
        target=target,
        controls=controls,
        rank=fit.rank,
        residual=residual,
        residual_q20=float(q20),
        residual_q80=float(q80),
        group=group,
        r_squared=fit.r_squared,
        coefficients=fit.coefficients,
    )


def _rank_control_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = pd.DataFrame(index=frame.index)
    for name in LOWER_CONTROLS_6:
        values = np.log(frame[name].to_numpy(dtype=float)) if name == "variance" else frame[name].to_numpy(dtype=float)
        result[f"rank_{name}"] = _unit_rank(values)
    return result


def _balance_rows(frame: pd.DataFrame, result: TargetResult, rank_controls: pd.DataFrame) -> pd.DataFrame:
    """Report raw and rank-scale balance without implying exact matching."""

    low = result.group == "low_q20"
    high = result.group == "high_q80"
    rows: list[dict[str, Any]] = []
    for name in result.controls:
        for representation, values in (
            ("raw", np.log(frame[name].to_numpy(float)) if name == "variance" else frame[name].to_numpy(float)),
            ("rank", rank_controls[f"rank_{name}"].to_numpy(float)),
        ):
            low_values = values[low]
            high_values = values[high]
            pooled = np.sqrt((np.var(low_values, ddof=1) + np.var(high_values, ddof=1)) / 2.0)
            smd = float((high_values.mean() - low_values.mean()) / pooled) if pooled > 0.0 else float("nan")
            rows.append(
                {
                    "dataset_id": result.dataset_id,
                    "target": result.target,
                    "control": _control_label(name),
                    "source_control": name,
                    "representation": representation,
                    "low_q20_mean": float(low_values.mean()),
                    "high_q80_mean": float(high_values.mean()),
                    "standardized_mean_difference": smd,
                    "low_q20_n": int(low.sum()),
                    "high_q80_n": int(high.sum()),
                }
            )
    return pd.DataFrame(rows)


def _ecdf_from_windows(standardized_values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if not len(standardized_values):
        return np.full(len(grid), np.nan, dtype=float)
    flattened = np.sort(standardized_values.reshape(-1))
    return np.searchsorted(flattened, grid, side="right") / float(len(flattened))


@dataclass(frozen=True)
class BlockAggregates:
    dataset_id: str
    target: str
    metric_names: tuple[str, ...]
    window_counts: np.ndarray
    metric_sums: np.ndarray
    ecdf_counts: np.ndarray


def _block_aggregates(
    prepared: PreparedDataset,
    result: TargetResult,
    metric_frame: pd.DataFrame,
    block_shape: tuple[int, int, int],
) -> BlockAggregates:
    """Cache group-specific sufficient statistics for block resampling."""

    block_labels = build_spatial_blocks(prepared.frame, block_shape)
    _, block_index = np.unique(block_labels, return_inverse=True)
    block_count = int(block_index.max()) + 1
    group_code = np.full(len(prepared.frame), -1, dtype=int)
    group_code[result.group == "low_q20"] = 0
    group_code[result.group == "high_q80"] = 1
    selected = group_code >= 0
    if not selected.any() or np.any(np.bincount(group_code[selected], minlength=2) == 0):
        raise ValueError(f"Residual groups are incomplete for {prepared.dataset_id}/{result.target}")
    metric_names = tuple(metric_frame.columns)
    window_counts = np.zeros((block_count, 2), dtype=np.int64)
    np.add.at(window_counts, (block_index[selected], group_code[selected]), 1)
    metric_sums = np.zeros((block_count, 2, len(metric_names)), dtype=float)
    for metric_index, name in enumerate(metric_names):
        np.add.at(metric_sums[:, :, metric_index], (block_index[selected], group_code[selected]), metric_frame.loc[selected, name].to_numpy(float))

    # Each uniform window contains exactly 75 values.  Accumulating a histogram
    # by spatial block makes a block-bootstrap ECDF inexpensive and preserves
    # equal window weights.
    sample_blocks = np.repeat(block_index[selected], prepared.standardized_values.shape[1])
    sample_groups = np.repeat(group_code[selected], prepared.standardized_values.shape[1])
    sample_values = prepared.standardized_values[selected].reshape(-1)
    first_grid_index = np.searchsorted(ECDF_GRID, sample_values, side="left")
    histogram = np.zeros((block_count, 2, len(ECDF_GRID) + 1), dtype=np.int64)
    np.add.at(histogram, (sample_blocks, sample_groups, first_grid_index), 1)
    ecdf_counts = np.cumsum(histogram[:, :, :-1], axis=2)
    return BlockAggregates(prepared.dataset_id, result.target, metric_names, window_counts, metric_sums, ecdf_counts)


def _max_abs_residual_control_correlation(residual: np.ndarray, design: np.ndarray) -> float:
    """Check that an OLS residual is numerically orthogonal to each control."""

    if design.shape[1] <= 1:
        return 0.0
    residual_sd = float(np.std(residual, ddof=0))
    if residual_sd <= 0.0 or not np.isfinite(residual_sd):
        raise ValueError("Rank residual has zero or non-finite variation")
    correlations: list[float] = []
    for index in range(1, design.shape[1]):
        control = design[:, index]
        control_sd = float(np.std(control, ddof=0))
        if control_sd <= 0.0 or not np.isfinite(control_sd):
            raise ValueError("Rank control has zero or non-finite variation")
        correlations.append(float(np.corrcoef(residual, control)[0, 1]))
    return float(np.max(np.abs(correlations))) if correlations else 0.0


def _association_block_aggregates(
    frame: pd.DataFrame,
    target_result: TargetResult,
    metric_frame: pd.DataFrame,
    rank_controls: pd.DataFrame,
    block_shape: tuple[int, int, int],
) -> tuple[AssociationBlockAggregates, pd.DataFrame, dict[str, RankOLSResult]]:
    """Prepare all rank-OLS association statistics by centre-coordinate block.

    Both the cumulant target and every empirical-distribution descriptor are
    ranked and projected on the same controls.  The stored matrices permit a
    block bootstrap to re-fit these projections in every replicate while
    keeping the full-subvolume rank transformation fixed.
    """

    controls = target_result.controls
    metric_names = tuple(str(name) for name in metric_frame.columns)
    if not metric_names:
        raise ValueError("Conditional association requires at least one descriptor")
    missing_primary = set(PRIMARY_ASSOCIATION_METRICS[target_result.target]).difference(metric_names)
    if missing_primary:
        raise ValueError(f"Tail metric frame omits primary descriptors: {sorted(missing_primary)}")
    design = _rank_design(rank_controls, controls)
    metric_fits = {
        name: _rank_response_projection(metric_frame[name].to_numpy(dtype=float), rank_controls, controls)
        for name in metric_names
    }
    z = np.column_stack((target_result.rank, *(metric_fits[name].rank for name in metric_names)))
    block_labels = build_spatial_blocks(frame, block_shape)
    _, block_index = np.unique(block_labels, return_inverse=True)
    block_count = int(block_index.max()) + 1
    design_width = design.shape[1]
    response_width = z.shape[1]
    xtx = np.zeros((block_count, design_width, design_width), dtype=float)
    xtz = np.zeros((block_count, design_width, response_width), dtype=float)
    ztz = np.zeros((block_count, response_width, response_width), dtype=float)
    for block_id in range(block_count):
        members = block_index == block_id
        block_design = design[members]
        block_response = z[members]
        xtx[block_id] = block_design.T @ block_design
        xtz[block_id] = block_design.T @ block_response
        ztz[block_id] = block_response.T @ block_response

    target_orthogonality = _max_abs_residual_control_correlation(target_result.residual, design)
    diagnostic_rows: list[dict[str, Any]] = []
    for name in metric_names:
        fit = metric_fits[name]
        diagnostic_rows.append(
            {
                "dataset_id": target_result.dataset_id,
                "target": target_result.target,
                "metric": name,
                "controls": ";".join(_control_label(control) for control in controls),
                "rank_transform_scope": "complete_subvolume_fixed",
                "target_rank_r_squared": float(target_result.r_squared),
                "metric_rank_r_squared": float(fit.r_squared),
                "target_max_abs_residual_control_correlation": target_orthogonality,
                "metric_max_abs_residual_control_correlation": _max_abs_residual_control_correlation(fit.residual, design),
                "centre_coordinate_block_count": int(block_count),
                "rank_design_width": int(design_width),
                "rank_design_condition_number": float(np.linalg.cond(design)),
            }
        )
    return (
        AssociationBlockAggregates(
            dataset_id=target_result.dataset_id,
            target=target_result.target,
            controls=controls,
            metric_names=metric_names,
            xtx=xtx,
            xtz=xtz,
            ztz=ztz,
        ),
        pd.DataFrame(diagnostic_rows),
        metric_fits,
    )


def _association_from_sufficient_statistics(
    xtx: np.ndarray,
    xtz: np.ndarray,
    ztz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return target-to-metric residual correlations from block sufficient statistics."""

    if xtx.ndim != 2 or xtx.shape[0] != xtx.shape[1]:
        raise ValueError("Rank-design cross-product is not square")
    if xtz.shape[0] != xtx.shape[0] or ztz.shape[0] != ztz.shape[1] or ztz.shape[0] != xtz.shape[1]:
        raise ValueError("Conditional-association sufficient statistics have incompatible shapes")
    if np.linalg.matrix_rank(xtx) != xtx.shape[0]:
        raise ValueError("Block-resampled rank design is singular")
    coefficients, _, _, _ = np.linalg.lstsq(xtx, xtz, rcond=None)
    residual_crossproducts = ztz - xtz.T @ coefficients
    residual_crossproducts = (residual_crossproducts + residual_crossproducts.T) / 2.0
    residual_ss = np.diag(residual_crossproducts)
    scale = max(float(np.max(np.abs(ztz))), 1.0)
    if not np.isfinite(residual_crossproducts).all() or np.any(residual_ss <= np.finfo(float).eps * scale):
        raise ValueError("Conditional rank residual has insufficient finite variation")
    target_ss = float(residual_ss[0])
    covariance = residual_crossproducts[0, 1:]
    metric_ss = residual_ss[1:]
    association = covariance / np.sqrt(target_ss * metric_ss)
    slope = covariance / target_ss
    if not np.isfinite(association).all() or not np.isfinite(slope).all() or np.any(np.abs(association) > 1.0 + 1.0e-10):
        raise ValueError("Conditional rank association is numerically invalid")
    return association, slope, residual_ss, float(np.linalg.cond(xtx))


def _point_association_rows(
    aggregates: list[AssociationBlockAggregates],
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Calculate per-subvolume and equal-subvolume point associations."""

    if not aggregates:
        raise ValueError("No association aggregates were supplied")
    target = aggregates[0].target
    metric_names = aggregates[0].metric_names
    controls = aggregates[0].controls
    if any(item.target != target or item.metric_names != metric_names or item.controls != controls for item in aggregates):
        raise ValueError("Association aggregate collection mixes incompatible definitions")
    rows: list[dict[str, Any]] = []
    point_values: dict[str, np.ndarray] = {}
    point_slopes: dict[str, np.ndarray] = {}
    for item in aggregates:
        association, slope, residual_ss, condition_number = _association_from_sufficient_statistics(
            item.xtx.sum(axis=0), item.xtz.sum(axis=0), item.ztz.sum(axis=0)
        )
        point_values[item.dataset_id] = association
        point_slopes[item.dataset_id] = slope
        for index, metric in enumerate(metric_names):
            rows.append(
                {
                    "target": target,
                    "metric": metric,
                    "dataset_id": item.dataset_id,
                    "controls": ";".join(_control_label(control) for control in controls),
                    "rank_transform_scope": "complete_subvolume_fixed",
                    "conditional_rank_association": float(association[index]),
                    "conditional_residual_slope": float(slope[index]),
                    "target_residual_rank_sum_squares": float(residual_ss[0]),
                    "metric_residual_rank_sum_squares": float(residual_ss[index + 1]),
                    "rank_design_condition_number": condition_number,
                }
            )
    point_values["macro_equal_dataset"] = np.mean(np.stack(list(point_values.values())), axis=0)
    point_slopes["macro_equal_dataset"] = np.mean(np.stack(list(point_slopes.values())), axis=0)
    for index, metric in enumerate(metric_names):
        rows.append(
            {
                "target": target,
                "metric": metric,
                "dataset_id": "macro_equal_dataset",
                "controls": ";".join(_control_label(control) for control in controls),
                "rank_transform_scope": "complete_subvolume_fixed",
                "conditional_rank_association": float(point_values["macro_equal_dataset"][index]),
                "conditional_residual_slope": float(point_slopes["macro_equal_dataset"][index]),
                "target_residual_rank_sum_squares": float("nan"),
                "metric_residual_rank_sum_squares": float("nan"),
                "rank_design_condition_number": float("nan"),
            }
        )
    return pd.DataFrame(rows), point_values


def _association_block_bootstrap(
    aggregates: list[AssociationBlockAggregates],
    *,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Bootstrap fixed-rank conditional associations by centre-coordinate blocks."""

    if replicates <= 0:
        raise ValueError("Bootstrap replicates must be positive")
    point_rows, point_values = _point_association_rows(aggregates)
    target = aggregates[0].target
    metric_names = aggregates[0].metric_names
    generator = np.random.default_rng(seed)
    draw_values: dict[str, list[np.ndarray]] = {item.dataset_id: [] for item in aggregates}
    draw_slopes: dict[str, list[np.ndarray]] = {item.dataset_id: [] for item in aggregates}
    macro_values: list[np.ndarray] = []
    macro_slopes: list[np.ndarray] = []
    accepted = 0
    attempts = 0
    max_attempts = replicates * 20
    while accepted < replicates:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError("Conditional-association block bootstrap could not obtain a full-rank replicate")
        replicate_values: list[np.ndarray] = []
        replicate_slopes: list[np.ndarray] = []
        current: list[tuple[str, np.ndarray, np.ndarray]] = []
        try:
            for item in aggregates:
                choices = generator.integers(0, len(item.xtx), size=len(item.xtx))
                association, slope, _, _ = _association_from_sufficient_statistics(
                    item.xtx[choices].sum(axis=0),
                    item.xtz[choices].sum(axis=0),
                    item.ztz[choices].sum(axis=0),
                )
                current.append((item.dataset_id, association, slope))
                replicate_values.append(association)
                replicate_slopes.append(slope)
        except ValueError:
            continue
        for dataset_id, association, slope in current:
            draw_values[dataset_id].append(association)
            draw_slopes[dataset_id].append(slope)
        macro_values.append(np.mean(np.stack(replicate_values), axis=0))
        macro_slopes.append(np.mean(np.stack(replicate_slopes), axis=0))
        accepted += 1

    sample_values = {dataset_id: np.stack(values) for dataset_id, values in draw_values.items()}
    sample_slopes = {dataset_id: np.stack(values) for dataset_id, values in draw_slopes.items()}
    sample_values["macro_equal_dataset"] = np.stack(macro_values)
    sample_slopes["macro_equal_dataset"] = np.stack(macro_slopes)
    summary_rows: list[dict[str, Any]] = []
    draw_rows: list[dict[str, Any]] = []
    for dataset_id, values in sample_values.items():
        slopes = sample_slopes[dataset_id]
        for metric_index, metric in enumerate(metric_names):
            summary_rows.append(
                {
                    "target": target,
                    "metric": metric,
                    "dataset_id": dataset_id,
                    "conditional_rank_association": float(point_values[dataset_id][metric_index]),
                    "ci95_low": float(np.quantile(values[:, metric_index], 0.025)),
                    "ci95_high": float(np.quantile(values[:, metric_index], 0.975)),
                    "conditional_residual_slope": float(
                        point_rows.loc[
                            (point_rows["dataset_id"] == dataset_id) & (point_rows["metric"] == metric),
                            "conditional_residual_slope",
                        ].iloc[0]
                    ),
                    "slope_ci95_low": float(np.quantile(slopes[:, metric_index], 0.025)),
                    "slope_ci95_high": float(np.quantile(slopes[:, metric_index], 0.975)),
                    "bootstrap_replicates": int(accepted),
                    "bootstrap_seed": int(seed),
                    "centre_coordinate_block_bootstrap": True,
                    "rank_transform_scope": "complete_subvolume_fixed",
                }
            )
        for replicate_index, (association_row, slope_row) in enumerate(zip(values, slopes, strict=True), start=1):
            for metric_index, metric in enumerate(metric_names):
                draw_rows.append(
                    {
                        "target": target,
                        "metric": metric,
                        "dataset_id": dataset_id,
                        "replicate": int(replicate_index),
                        "conditional_rank_association": float(association_row[metric_index]),
                        "conditional_residual_slope": float(slope_row[metric_index]),
                        "bootstrap_seed": int(seed),
                    }
                )
    return pd.DataFrame(summary_rows), pd.DataFrame(draw_rows), {"accepted": accepted, "attempts": attempts}


def _point_summary(aggregates: BlockAggregates) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = aggregates.window_counts.sum(axis=0)
    if np.any(counts == 0):
        raise ValueError(f"Point support has an absent residual group for {aggregates.dataset_id}/{aggregates.target}")
    metrics = aggregates.metric_sums.sum(axis=0) / counts[:, None]
    ecdf = aggregates.ecdf_counts.sum(axis=0) / (counts[:, None] * WINDOW[0] * WINDOW[1] * WINDOW[2])
    return counts, metrics, ecdf


def _block_bootstrap(
    aggregates: list[BlockAggregates],
    *,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Return block-bootstrap CIs for per-volume and equal-volume contrasts."""

    if replicates <= 0:
        raise ValueError("Bootstrap replicates must be positive")
    target = aggregates[0].target
    if any(item.target != target for item in aggregates):
        raise ValueError("Bootstrap collection mixes cumulant targets")
    metric_names = aggregates[0].metric_names
    if any(item.metric_names != metric_names for item in aggregates):
        raise ValueError("Bootstrap collection mixes metric definitions")
    point: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {
        item.dataset_id: _point_summary(item) for item in aggregates
    }
    generator = np.random.default_rng(seed)
    metric_samples: dict[str, list[np.ndarray]] = {item.dataset_id: [] for item in aggregates}
    ecdf_samples: dict[str, list[np.ndarray]] = {item.dataset_id: [] for item in aggregates}
    macro_metrics: list[np.ndarray] = []
    macro_ecdf: list[np.ndarray] = []
    accepted = 0
    attempts = 0
    max_attempts = replicates * 100
    while accepted < replicates:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError("Block bootstrap could not obtain complete residual-group support")
        replicate_metrics: list[np.ndarray] = []
        replicate_ecdf: list[np.ndarray] = []
        per_dataset: list[tuple[str, np.ndarray, np.ndarray]] = []
        valid = True
        for item in aggregates:
            choices = generator.integers(0, len(item.window_counts), size=len(item.window_counts))
            counts = item.window_counts[choices].sum(axis=0)
            if np.any(counts == 0):
                valid = False
                break
            means = item.metric_sums[choices].sum(axis=0) / counts[:, None]
            ecdf = item.ecdf_counts[choices].sum(axis=0) / (counts[:, None] * WINDOW[0] * WINDOW[1] * WINDOW[2])
            contrast_metrics = means[1] - means[0]
            contrast_ecdf = ecdf[1] - ecdf[0]
            per_dataset.append((item.dataset_id, contrast_metrics, contrast_ecdf))
            replicate_metrics.append(contrast_metrics)
            replicate_ecdf.append(contrast_ecdf)
        if not valid:
            continue
        for dataset_id, contrast_metrics, contrast_ecdf in per_dataset:
            metric_samples[dataset_id].append(contrast_metrics)
            ecdf_samples[dataset_id].append(contrast_ecdf)
        macro_metrics.append(np.mean(np.stack(replicate_metrics), axis=0))
        macro_ecdf.append(np.mean(np.stack(replicate_ecdf), axis=0))
        accepted += 1

    metric_rows: list[dict[str, Any]] = []
    ecdf_rows: list[dict[str, Any]] = []
    point_metric_contrasts = {identifier: values[1][1] - values[1][0] for identifier, values in point.items()}
    point_ecdf_contrasts = {identifier: values[2][1] - values[2][0] for identifier, values in point.items()}
    sample_sets = {
        **{identifier: (np.stack(metric_samples[identifier]), np.stack(ecdf_samples[identifier])) for identifier in metric_samples},
        "macro_equal_dataset": (np.stack(macro_metrics), np.stack(macro_ecdf)),
    }
    point_metric_contrasts["macro_equal_dataset"] = np.mean(np.stack(list(point_metric_contrasts.values())), axis=0)
    point_ecdf_contrasts["macro_equal_dataset"] = np.mean(np.stack(list(point_ecdf_contrasts.values())), axis=0)
    for dataset_id, (metric_draws, ecdf_draws) in sample_sets.items():
        for index, name in enumerate(metric_names):
            metric_rows.append(
                {
                    "target": target,
                    "dataset_id": dataset_id,
                    "metric": name,
                    "contrast_high_q80_minus_low_q20": float(point_metric_contrasts[dataset_id][index]),
                    "ci95_low": float(np.quantile(metric_draws[:, index], 0.025)),
                    "ci95_high": float(np.quantile(metric_draws[:, index], 0.975)),
                    "bootstrap_replicates": int(accepted),
                    "bootstrap_seed": int(seed),
                }
            )
        for index, threshold in enumerate(ECDF_GRID):
            ecdf_rows.append(
                {
                    "target": target,
                    "dataset_id": dataset_id,
                    "standardized_porosity": float(threshold),
                    "ecdf_difference_high_q80_minus_low_q20": float(point_ecdf_contrasts[dataset_id][index]),
                    "ci95_low": float(np.quantile(ecdf_draws[:, index], 0.025)),
                    "ci95_high": float(np.quantile(ecdf_draws[:, index], 0.975)),
                    "bootstrap_replicates": int(accepted),
                    "bootstrap_seed": int(seed),
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(ecdf_rows), {"accepted": accepted, "attempts": attempts}


def _point_rows(aggregates: list[BlockAggregates]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, Any]] = []
    ecdf_rows: list[dict[str, Any]] = []
    per_dataset_metrics: list[np.ndarray] = []
    per_dataset_ecdf: list[np.ndarray] = []
    for item in aggregates:
        counts, metrics, ecdf = _point_summary(item)
        per_dataset_metrics.append(metrics)
        per_dataset_ecdf.append(ecdf)
        for group_index, group in enumerate(("low_q20", "high_q80")):
            for metric_index, name in enumerate(item.metric_names):
                metric_rows.append(
                    {
                        "target": item.target,
                        "dataset_id": item.dataset_id,
                        "group": group,
                        "metric": name,
                        "value": float(metrics[group_index, metric_index]),
                        "window_count": int(counts[group_index]),
                    }
                )
            for threshold, value in zip(ECDF_GRID, ecdf[group_index], strict=True):
                ecdf_rows.append(
                    {
                        "target": item.target,
                        "dataset_id": item.dataset_id,
                        "group": group,
                        "standardized_porosity": float(threshold),
                        "ecdf": float(value),
                        "window_count": int(counts[group_index]),
                    }
                )
    macro_metrics = np.mean(np.stack(per_dataset_metrics), axis=0)
    macro_ecdf = np.mean(np.stack(per_dataset_ecdf), axis=0)
    target = aggregates[0].target
    names = aggregates[0].metric_names
    for group_index, group in enumerate(("low_q20", "high_q80")):
        for metric_index, name in enumerate(names):
            metric_rows.append(
                {
                    "target": target,
                    "dataset_id": "macro_equal_dataset",
                    "group": group,
                    "metric": name,
                    "value": float(macro_metrics[group_index, metric_index]),
                    "window_count": int(sum(_point_summary(item)[0][group_index] for item in aggregates)),
                }
            )
        for threshold, value in zip(ECDF_GRID, macro_ecdf[group_index], strict=True):
            ecdf_rows.append(
                {
                    "target": target,
                    "dataset_id": "macro_equal_dataset",
                    "group": group,
                    "standardized_porosity": float(threshold),
                    "ecdf": float(value),
                    "window_count": int(sum(_point_summary(item)[0][group_index] for item in aggregates)),
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(ecdf_rows)


def _windows_disjoint(high: pd.Series, low: pd.Series) -> bool:
    """Return whether two centered 5x5x3 supports share no source voxel."""

    return (
        abs(int(high.x) - int(low.x)) >= WINDOW[0]
        or abs(int(high.y) - int(low.y)) >= WINDOW[1]
        or abs(int(high.z) - int(low.z)) >= WINDOW[2]
    )


def illustrative_pair(
    frame: pd.DataFrame,
    result: TargetResult,
    rank_controls: pd.DataFrame,
    *,
    caliper_iqr: float,
) -> tuple[dict[str, Any], tuple[int, int] | None]:
    """Select one deterministic, non-overlapping illustrative conditional pair.

    The pair is only for the window-level ECDF illustration.  Quantitative
    inference comes from the full residual groups and spatial-block bootstrap.
    """

    low_index = np.flatnonzero(result.group == "low_q20")
    high_index = np.flatnonzero(result.group == "high_q80")
    controls = list(result.controls)
    metadata: dict[str, Any] = {
        "dataset_id": result.dataset_id,
        "target": result.target,
        "controls": ";".join(controls),
        "caliper_iqr": float(caliper_iqr),
        "high_anchor_count": int(len(high_index)),
        "low_candidate_count": int(len(low_index)),
        "nearest_candidates_examined_per_high_anchor": int(min(ILLUSTRATIVE_PAIR_NEIGHBORS, len(low_index))),
    }
    if not len(low_index) or not len(high_index):
        metadata.update(
            {
                "high_anchor_coverage_within_examined_neighbours": float("nan"),
                "eligible_nonoverlapping_pair_count": 0,
                "selected": False,
            }
        )
        return metadata, None
    controls_rank = rank_controls.loc[:, [f"rank_{name}" for name in controls]].to_numpy(float)
    iqr = np.quantile(controls_rank, 0.75, axis=0) - np.quantile(controls_rank, 0.25, axis=0)
    iqr = np.where(iqr > 0.0, iqr, 1.0)
    low_points = controls_rank[low_index] / iqr
    high_points = controls_rank[high_index] / iqr
    tree = cKDTree(low_points)
    # These pairs are traceable illustrations only. Restricting to a fixed
    # number of nearest candidates prevents pair drawing from dominating the
    # all-window ECDF/bootstrap computation that provides the actual evidence.
    neighbours = min(ILLUSTRATIVE_PAIR_NEIGHBORS, len(low_index))
    distances, local_indices = tree.query(high_points, k=neighbours)
    if neighbours == 1:
        distances = distances[:, None]
        local_indices = local_indices[:, None]
    high_has_candidate = np.zeros(len(high_index), dtype=bool)
    options: list[tuple[float, float, int, int]] = []
    for high_position, (distance_row, low_local_row) in enumerate(zip(distances, local_indices, strict=True)):
        high_row = int(high_index[high_position])
        for distance, low_local in zip(np.ravel(distance_row), np.ravel(low_local_row), strict=True):
            low_row = int(low_index[int(low_local)])
            difference = np.abs(controls_rank[high_row] - controls_rank[low_row]) / iqr
            chebyshev = float(difference.max())
            if chebyshev > caliper_iqr:
                continue
            high_has_candidate[high_position] = True
            if _windows_disjoint(frame.iloc[high_row], frame.iloc[low_row]):
                options.append((chebyshev, float(distance), high_row, low_row))
    coverage = float(high_has_candidate.mean()) if len(high_has_candidate) else float("nan")
    metadata.update(
        {
            "high_anchor_coverage_within_examined_neighbours": coverage,
            "eligible_nonoverlapping_pair_count": int(len(options)),
        }
    )
    if not options:
        metadata["selected"] = False
        return metadata, None
    # A fixed ordering makes the selection reproducible and avoids visual
    # selection based on tail summaries or any external label.
    _, _, high_row, low_row = min(options, key=lambda item: (item[0], item[1], item[2], item[3]))
    high = frame.iloc[high_row]
    low = frame.iloc[low_row]
    metadata.update(
        {
            "selected": True,
            "high_row": int(high_row),
            "low_row": int(low_row),
            "high_x": int(high.x),
            "high_y": int(high.y),
            "high_z": int(high.z),
            "low_x": int(low.x),
            "low_y": int(low.y),
            "low_z": int(low.z),
            "high_target": float(high[result.target]),
            "low_target": float(low[result.target]),
            "high_residual": float(result.residual[high_row]),
            "low_residual": float(result.residual[low_row]),
        }
    )
    for name in controls:
        high_value = float(high[name])
        low_value = float(low[name])
        metadata[f"high_{name}"] = high_value
        metadata[f"low_{name}"] = low_value
        metadata[f"difference_{name}"] = high_value - low_value
    return metadata, (high_row, low_row)


def _plot_ecdf(
    ecdf_points: pd.DataFrame,
    ecdf_contrasts: pd.DataFrame,
    target: str,
    output_path: Path,
) -> None:
    datasets = ["z15_26", "z45_56", "z75_86"]
    figure, axes = plt.subplots(2, len(datasets), figsize=(13.2, 5.8), sharex=True, constrained_layout=True)
    for column, dataset_id in enumerate(datasets):
        current = ecdf_points.loc[(ecdf_points["target"] == target) & (ecdf_points["dataset_id"] == dataset_id)]
        upper = axes[0, column]
        lower = axes[1, column]
        for group, color in (("low_q20", "#277da1"), ("high_q80", "#d1495b")):
            data = current.loc[current["group"] == group].sort_values("standardized_porosity")
            upper.step(data["standardized_porosity"], data["ecdf"], where="post", color=color, label=group)
        upper.set_title(dataset_id)
        upper.set_ylabel("mean local ECDF")
        upper.grid(alpha=0.25)
        difference = ecdf_contrasts.loc[
            (ecdf_contrasts["target"] == target) & (ecdf_contrasts["dataset_id"] == dataset_id)
        ].sort_values("standardized_porosity")
        lower.plot(difference["standardized_porosity"], difference["ecdf_difference_high_q80_minus_low_q20"], color="#5d3a9b")
        lower.fill_between(
            difference["standardized_porosity"], difference["ci95_low"], difference["ci95_high"],
            color="#5d3a9b", alpha=0.22,
        )
        lower.axhline(0.0, color="black", linewidth=0.8)
        lower.set(xlabel="within-window standardized porosity", ylabel="high - low ECDF")
        lower.grid(alpha=0.25)
    axes[0, 0].legend(frameon=False, loc="best")
    figure.suptitle(f"Conditional local empirical distributions: {target}")
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def _plot_selected_pairs(
    pair_rows: pd.DataFrame,
    prepared: Mapping[str, PreparedDataset],
    selected_indices: Mapping[tuple[str, str], tuple[int, int]],
    output_path: Path,
) -> None:
    datasets = ["z15_26", "z45_56", "z75_86"]
    targets = ["lambda5", "lambda6"]
    figure, axes = plt.subplots(2, 3, figsize=(13.2, 6.6), sharex=True, sharey=True, constrained_layout=True)
    for row_index, target in enumerate(targets):
        for column, dataset_id in enumerate(datasets):
            axis = axes[row_index, column]
            selected = selected_indices.get((dataset_id, target))
            if selected is None:
                axis.text(0.5, 0.5, "No eligible pair", ha="center", va="center", transform=axis.transAxes)
                continue
            high_index, low_index = selected
            values = prepared[dataset_id].standardized_values
            low = np.sort(values[low_index])
            high = np.sort(values[high_index])
            low_y = np.arange(1, len(low) + 1) / len(low)
            high_y = np.arange(1, len(high) + 1) / len(high)
            axis.step(low, low_y, where="post", color="#277da1", label="low residual")
            axis.step(high, high_y, where="post", color="#d1495b", label="high residual")
            axis.set_title(f"{dataset_id}, {target}")
            axis.grid(alpha=0.25)
            if row_index == len(targets) - 1:
                axis.set_xlabel("within-window standardized porosity")
            if column == 0:
                axis.set_ylabel("window ECDF")
    axes[0, 0].legend(frameon=False, loc="best")
    figure.suptitle("Deterministic illustrative residual-window pairs")
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def _plot_tail_contrasts(contrasts: pd.DataFrame, output_path: Path) -> None:
    selected = [
        ("lambda5", "tail_imbalance", "lambda5: tail imbalance"),
        ("lambda6", "outer_tail_probability", "lambda6: outer-tail probability"),
        ("lambda6", "outer_span_q05_q95", "lambda6: q95 - q05"),
    ]
    figure, axes = plt.subplots(1, len(selected), figsize=(12.5, 3.7), constrained_layout=True)
    macro = contrasts.loc[contrasts["dataset_id"] == "macro_equal_dataset"]
    for axis, (target, metric, title) in zip(axes, selected, strict=True):
        row = macro.loc[(macro["target"] == target) & (macro["metric"] == metric)]
        if row.empty:
            axis.text(0.5, 0.5, "Unavailable", ha="center", va="center", transform=axis.transAxes)
            continue
        value = row.iloc[0]
        estimate = float(value["contrast_high_q80_minus_low_q20"])
        low = float(value["ci95_low"])
        high = float(value["ci95_high"])
        axis.errorbar([0], [estimate], yerr=[[estimate - low], [high - estimate]], fmt="o", color="#5d3a9b", capsize=4)
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set(xlim=(-0.8, 0.8), xticks=[], title=title, ylabel="high residual - low residual")
        axis.grid(axis="y", alpha=0.25)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def _plot_conditional_rank_associations(associations: pd.DataFrame, output_path: Path) -> None:
    """Render the pre-specified conditional associations as a compact forest plot."""

    selected = [
        ("lambda5", "tail_imbalance", "lambda5 and tail imbalance"),
        ("lambda6", "outer_tail_probability", "lambda6 and outer-tail probability"),
        ("lambda6", "outer_span_q05_q95", "lambda6 and Q95 - Q05"),
    ]
    dataset_order = ["z15_26", "z45_56", "z75_86", "macro_equal_dataset"]
    label_map = {"z15_26": "z15-26", "z45_56": "z45-56", "z75_86": "z75-86", "macro_equal_dataset": "equal-volume mean"}
    figure, axes = plt.subplots(1, len(selected), figsize=(12.8, 4.2), sharey=True, constrained_layout=True)
    for axis, (target, metric, title) in zip(axes, selected, strict=True):
        rows = associations.loc[
            (associations["target"] == target) & (associations["metric"] == metric)
        ].set_index("dataset_id").reindex(dataset_order)
        if rows["conditional_rank_association"].isna().any():
            raise ValueError(f"Conditional association figure lacks rows for {target}/{metric}")
        y = np.arange(len(dataset_order), dtype=float)
        estimates = rows["conditional_rank_association"].to_numpy(dtype=float)
        low = rows["ci95_low"].to_numpy(dtype=float)
        high = rows["ci95_high"].to_numpy(dtype=float)
        colors = ["#277da1", "#277da1", "#277da1", "#d1495b"]
        for position, estimate, lower, upper, color in zip(y, estimates, low, high, colors, strict=True):
            axis.errorbar(
                estimate,
                position,
                xerr=[[estimate - lower], [upper - estimate]],
                fmt="o",
                color=color,
                capsize=3,
            )
        finite_bounds = np.concatenate((low, high, np.asarray([0.0])))
        extent = max(float(np.max(np.abs(finite_bounds))), 0.05)
        axis.axvline(0.0, color="black", linewidth=0.8)
        axis.set(
            xlim=(-1.1 * extent, 1.1 * extent),
            yticks=y,
            yticklabels=[label_map[item] for item in dataset_order],
            title=title,
            xlabel="fixed-rank conditional association",
        )
        axis.grid(axis="x", alpha=0.25)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def _format_table(frame: pd.DataFrame, columns: list[str], limit: int | None = None) -> str:
    if frame.empty:
        return "No rows."
    selected = frame.reindex(columns=columns).copy()
    if limit is not None:
        selected = selected.head(limit)
    # Avoid a report-only dependency on the optional ``tabulate`` package.
    def cell(value: Any) -> str:
        if value is None or pd.isna(value):
            return ""
        if isinstance(value, (bool, np.bool_)):
            return "true" if bool(value) else "false"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.4g}"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = ["| " + " | ".join(cell(value) for value in row) + " |" for row in selected.itertuples(index=False, name=None)]
    return "\n".join([header, separator, *rows])


def _write_report(
    path: Path,
    *,
    source_hashes: Mapping[str, Any],
    support: pd.DataFrame,
    residual_summary: pd.DataFrame,
    association_summary: pd.DataFrame,
    association_diagnostics: pd.DataFrame,
    contrasts: pd.DataFrame,
    pair_rows: pd.DataFrame,
    bootstrap_diagnostics: Mapping[str, Any],
) -> None:
    primary_associations = association_summary.loc[
        association_summary.apply(
            lambda row: str(row["metric"]) in PRIMARY_ASSOCIATION_METRICS.get(str(row["target"]), ()), axis=1
        )
        & association_summary["dataset_id"].eq("macro_equal_dataset")
    ].copy()
    primary_diagnostics = association_diagnostics.loc[
        association_diagnostics.apply(
            lambda row: str(row["metric"]) in PRIMARY_ASSOCIATION_METRICS.get(str(row["target"]), ()), axis=1
        )
    ].copy()
    macro_metrics = contrasts.loc[
        (contrasts["dataset_id"] == "macro_equal_dataset")
        & (
            ((contrasts["target"] == "lambda5") & (contrasts["metric"] == "tail_imbalance"))
            | ((contrasts["target"] == "lambda6") & contrasts["metric"].isin(["outer_tail_probability", "outer_span_q05_q95"]))
        )
    ].copy()
    lines = [
        "# Conditional High-Order Window Analysis",
        "",
        "## Scope",
        "",
        "This post-processing analysis uses only the frozen run-local porosity grids and uniform 5x5x3 Section 4.5 numerical maps. Channel labels, facies validity, prediction models, and C++ reruns are not used.",
        "",
        "For `lambda5`, the rank projection conditions on mean, log variance, high-porosity proportion, lambda3, and lambda4. For `lambda6`, it additionally conditions on lambda5. The Q20 and Q80 groups are defined from the resulting rank-OLS residuals within each subvolume.",
        "",
        "## Numerical Support",
        "",
        _format_table(
            support,
            [
                "dataset_id", "structural_rows", "numerical_windows", "window_samples",
                "reference_variance", "variance_epsilon", "minimum_variance_retained",
                "minimum_ess", "maximum_ess", "formal_analysis_variance_flag_disagreement",
            ],
        ),
        "",
        "Retained windows are required independently to have finite numerical fields, `cpp_valid=true`, effective sample size 75, and `kappa2 > variance_epsilon`. The accompanying `numerical_support.csv` records the exclusion counts for each gate.",
        "",
        "## Conditional Rank Variation",
        "",
        _format_table(
            residual_summary,
            ["dataset_id", "target", "controls", "r_squared", "unexplained_rank_variance", "low_q20_count", "high_q80_count"],
        ),
        "",
        "The unexplained rank variance is a descriptive non-redundancy measure. It is not a prediction score or an external geological-validation metric.",
        "",
        "## Fixed-Rank Conditional Associations",
        "",
        "For each empirical-distribution descriptor, its complete-subvolume rank and the cumulant rank are separately projected on the same lower-order rank controls. The reported association is the correlation of those two residuals. Within every centre-coordinate spatial-block bootstrap replicate, both projections are re-fit from the resampled block sufficient statistics; the full-subvolume rank transformation remains fixed.",
        "",
        _format_table(
            primary_associations,
            [
                "target", "metric", "conditional_rank_association", "ci95_low", "ci95_high",
                "conditional_residual_slope", "bootstrap_replicates",
            ],
        ),
        "",
        "The three rows above were specified before this analysis: lambda5 with tail imbalance, and lambda6 with outer-tail probability and Q95-Q05. `conditional_rank_associations.csv` contains the complete descriptor matrix. The local empirical-distribution descriptors and cumulants are calculated from the same windows, so this is a conditional shape-association analysis rather than an external prediction experiment.",
        "",
        _format_table(
            primary_diagnostics,
            [
                "dataset_id", "target", "metric", "target_rank_r_squared", "metric_rank_r_squared",
                "target_max_abs_residual_control_correlation", "metric_max_abs_residual_control_correlation",
                "centre_coordinate_block_count",
            ],
        ),
        "",
        "## Block-Bootstrap Local-Distribution Contrasts",
        "",
        _format_table(
            macro_metrics,
            ["target", "metric", "contrast_high_q80_minus_low_q20", "ci95_low", "ci95_high", "bootstrap_replicates"],
        ),
        "",
        "Each contrast is an equal-subvolume average. Spatial blocks, rather than individual samples, were resampled by window centre coordinate. The ECDF figures and Q20/Q80 contrasts remain descriptive group summaries; the fixed-rank conditional associations above provide the lower-order-adjusted quantitative comparison. Outer-tail summaries are compact descriptors rather than a claim that a single cumulant represents only one tail property.",
        "",
        "## Illustrative Window Pairs",
        "",
        _format_table(
            pair_rows,
            [
                "dataset_id", "target", "selected", "nearest_candidates_examined_per_high_anchor",
                "high_anchor_coverage_within_examined_neighbours", "eligible_nonoverlapping_pair_count",
                "high_x", "high_y", "high_z", "low_x", "low_y", "low_z",
            ],
        ),
        "",
        "The selected pairs are traceable illustrations. The quantitative evidence is based on all residual-group windows and the spatial-block bootstrap, not on a selected pair.",
        "",
        "## Bootstrap Diagnostics",
        "",
        "```json",
        json.dumps(bootstrap_diagnostics, indent=2, sort_keys=True),
        "```",
        "",
        "## Source Hashes",
        "",
        "```json",
        json.dumps(source_hashes, indent=2, sort_keys=True),
        "```",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_hash_manifest(output_directory: Path) -> None:
    manifest = output_directory / "sha256_manifest.csv"
    paths = sorted(
        (path for path in output_directory.rglob("*") if path.is_file() and path != manifest),
        key=lambda path: str(path.relative_to(output_directory)).lower(),
    )
    rows = [{"path": str(path.relative_to(output_directory)).replace("\\", "/"), "sha256": sha256(path), "bytes": path.stat().st_size} for path in paths]
    pd.DataFrame(rows).to_csv(manifest, index=False)


def run_analysis(arguments: argparse.Namespace) -> dict[str, Any]:
    run_directory = arguments.run_dir.resolve()
    output_directory = arguments.output_dir.resolve()
    if not (run_directory / "direct_analysis_values.csv").is_file():
        raise FileNotFoundError(f"Missing direct_analysis_values.csv under {run_directory}")
    if output_directory == run_directory or output_directory.is_relative_to(run_directory):
        raise ValueError("Output directory must be outside the completed formal run directory")
    if output_directory.exists() and any(output_directory.iterdir()):
        raise FileExistsError(f"Output directory must be new or empty: {output_directory}")
    formal_bundle = _verify_formal_run(run_directory)
    frozen = _load_frozen_inputs(run_directory, arguments.common_config, formal_bundle)
    block_shape = _formal_block_shape(frozen.common, arguments.block_shape)
    output_directory.mkdir(parents=True, exist_ok=True)
    figures_directory = output_directory / "figures"
    figures_directory.mkdir(exist_ok=True)

    values_path = Path(formal_bundle["direct_analysis_values_path"])
    selected, support_frame = _select_numeric_support(
        values_path,
        frozen.variance_epsilon,
        frozen.reference_variance,
    )
    source_hashes: dict[str, Any] = {
        "analysis_utc": datetime.now(timezone.utc).isoformat(),
        "formal_run_directory": str(run_directory),
        "direct_analysis_values": {
            "run_relative_path": str(values_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": str(formal_bundle["direct_analysis_values_sha256"]),
        },
        "run_manifest": {
            "run_relative_path": str(Path(formal_bundle["run_manifest_path"]).relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.run_manifest_sha256,
        },
        "formal_gate": {
            "run_relative_path": str(frozen.gate_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.gate_sha256,
        },
        "formal_evidence_manifest": {
            "run_relative_path": str(frozen.evidence_manifest_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.evidence_manifest_sha256,
        },
        "formal_sha256_inventory": {
            "run_relative_path": str(frozen.archived_hash_manifest_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.archived_hash_manifest_sha256,
        },
        "frozen_common_protocol": {
            "run_relative_path": str(frozen.common_config_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.common_config_sha256,
        },
        "frozen_section45_protocol": {
            "run_relative_path": str(Path(formal_bundle["direct_section45_path"]).relative_to(run_directory)).replace("\\", "/"),
            "sha256": str(formal_bundle["direct_section45_sha256"]),
        },
        "preprocessing_manifest": {
            "run_relative_path": str(frozen.preprocessing_manifest_path.relative_to(run_directory)).replace("\\", "/"),
            "sha256": frozen.preprocessing_manifest_sha256,
        },
        "run_local_porosity_grids": {
            identifier: {
                "run_relative_path": str(path.relative_to(run_directory)).replace("\\", "/"),
                "sha256": frozen.porosity_hashes[identifier],
                "reference_variance": float(frozen.reference_variance[identifier]),
                "variance_epsilon": float(frozen.variance_epsilon[identifier]),
            }
            for identifier, path in sorted(frozen.porosity_paths.items())
        },
        "high_porosity_threshold": float(frozen.common["high_porosity_threshold"]),
        "uniform_window": list(WINDOW),
        "required_effective_sample_size": UNIFORM_WINDOW_SAMPLE_SIZE,
        "spatial_block_shape_xyz": list(block_shape),
        "conditional_association_protocol": {
            "rank_transform_scope": "complete_subvolume_fixed",
            "bootstrap_refits_rank_ols": True,
            "bootstrap_resampling_unit": "centre_coordinate_spatial_block",
            "primary_metrics": {key: list(value) for key, value in PRIMARY_ASSOCIATION_METRICS.items()},
        },
        "script": str(SCRIPT_PATH),
        "script_sha256": sha256(SCRIPT_PATH),
    }

    prepared: dict[str, PreparedDataset] = {}
    parity_reports: list[pd.DataFrame] = []
    high_threshold = float(frozen.common["high_porosity_threshold"])
    for dataset_id, local_frame in selected.groupby("dataset_id", sort=True):
        ordered = local_frame.sort_values(["z", "y", "x"], kind="stable").reset_index(drop=True)
        epsilon = float(frozen.variance_epsilon[str(dataset_id)])
        values, _ = _recover_dataset_windows(ordered, frozen.porosity_paths[str(dataset_id)])
        standardized = _standardize_windows(values, epsilon)
        item = PreparedDataset(str(dataset_id), ordered, values, standardized)
        parity_reports.append(
            parity_audit(
                item,
                high_threshold=high_threshold,
                variance_epsilon=epsilon,
                sample_size=int(arguments.parity_sample_size),
                rtol=float(arguments.parity_rtol),
                atol=float(arguments.parity_atol),
            )
        )
        prepared[item.dataset_id] = item
    parity_frame = pd.concat(parity_reports, ignore_index=True)
    parity_frame.to_csv(output_directory / "reconstruction_parity_audit.csv", index=False)
    support_frame.to_csv(output_directory / "numerical_support.csv", index=False)

    residual_rows: list[dict[str, Any]] = []
    balance_frames: list[pd.DataFrame] = []
    value_frames: list[pd.DataFrame] = []
    pair_rows: list[dict[str, Any]] = []
    selected_pair_indices: dict[tuple[str, str], tuple[int, int]] = {}
    target_aggregates: dict[str, list[BlockAggregates]] = {"lambda5": [], "lambda6": []}
    target_association_aggregates: dict[str, list[AssociationBlockAggregates]] = {"lambda5": [], "lambda6": []}
    association_diagnostic_frames: list[pd.DataFrame] = []
    for dataset_id in sorted(prepared):
        item = prepared[dataset_id]
        local = item.frame.copy()
        rank_controls = _rank_control_frame(local)
        metrics = _tail_metrics(item.standardized_values)
        for target, controls in (("lambda5", LOWER_CONTROLS_5), ("lambda6", LOWER_CONTROLS_6)):
            result = _rank_projection(local, target, controls, rank_controls=rank_controls)
            local[f"{target}_rank_residual"] = result.residual
            local[f"{target}_residual_group"] = result.group
            for name in LOWER_CONTROLS_6:
                local[f"rank_{name}"] = rank_controls[f"rank_{name}"].to_numpy(float)
            residual_rows.append(
                {
                    "dataset_id": dataset_id,
                    "target": target,
                    "controls": ";".join(_control_label(control) for control in controls),
                    "r_squared": result.r_squared,
                    "unexplained_rank_variance": 1.0 - result.r_squared,
                    "residual_q20": result.residual_q20,
                    "residual_q80": result.residual_q80,
                    "low_q20_count": int(np.sum(result.group == "low_q20")),
                    "high_q80_count": int(np.sum(result.group == "high_q80")),
                    **{f"coefficient_{index}": float(value) for index, value in enumerate(result.coefficients)},
                }
            )
            balance_frames.append(_balance_rows(local, result, rank_controls))
            aggregate = _block_aggregates(
                PreparedDataset(dataset_id, local, item.values, item.standardized_values),
                result,
                metrics,
                block_shape,
            )
            target_aggregates[target].append(aggregate)
            association_aggregate, association_diagnostics, metric_fits = _association_block_aggregates(
                local,
                result,
                metrics,
                rank_controls,
                block_shape,
            )
            target_association_aggregates[target].append(association_aggregate)
            association_diagnostic_frames.append(association_diagnostics)
            for metric in PRIMARY_ASSOCIATION_METRICS[target]:
                local[f"{target}_{metric}_rank_residual"] = metric_fits[metric].residual
            pair_row, indices = illustrative_pair(local, result, rank_controls, caliper_iqr=float(arguments.pair_caliper_iqr))
            pair_rows.append(pair_row)
            if indices is not None:
                selected_pair_indices[(dataset_id, target)] = indices
        value_frames.append(local)

    residual_frame = pd.DataFrame(residual_rows).sort_values(["target", "dataset_id"], kind="stable")
    residual_frame.to_csv(output_directory / "conditional_residual_summary.csv", index=False)
    balance_frame = pd.concat(balance_frames, ignore_index=True)
    balance_frame.to_csv(output_directory / "conditional_residual_balance.csv", index=False)
    values_frame = pd.concat(value_frames, ignore_index=True)
    values_frame.to_csv(output_directory / "conditional_window_residuals.csv", index=False)
    pairs_frame = pd.DataFrame(pair_rows).sort_values(["target", "dataset_id"], kind="stable")
    pairs_frame.to_csv(output_directory / "illustrative_conditional_pairs.csv", index=False)

    point_metric_frames: list[pd.DataFrame] = []
    point_ecdf_frames: list[pd.DataFrame] = []
    contrast_frames: list[pd.DataFrame] = []
    ecdf_contrast_frames: list[pd.DataFrame] = []
    bootstrap_diagnostics: dict[str, Mapping[str, int]] = {}
    for target, aggregates in target_aggregates.items():
        point_metrics, point_ecdf = _point_rows(aggregates)
        contrasts, ecdf_contrasts, diagnostics = _block_bootstrap(
            aggregates,
            replicates=int(arguments.bootstrap_replicates),
            seed=int(arguments.seed) + (5 if target == "lambda5" else 6),
        )
        point_metric_frames.append(point_metrics)
        point_ecdf_frames.append(point_ecdf)
        contrast_frames.append(contrasts)
        ecdf_contrast_frames.append(ecdf_contrasts)
        bootstrap_diagnostics[target] = diagnostics
    point_metrics = pd.concat(point_metric_frames, ignore_index=True)
    point_ecdf = pd.concat(point_ecdf_frames, ignore_index=True)
    contrasts = pd.concat(contrast_frames, ignore_index=True)
    ecdf_contrasts = pd.concat(ecdf_contrast_frames, ignore_index=True)
    point_metrics.to_csv(output_directory / "conditional_tail_summary.csv", index=False)
    contrasts.to_csv(output_directory / "conditional_tail_contrasts_bootstrap.csv", index=False)
    point_ecdf.to_csv(output_directory / "conditional_ecdf.csv", index=False)
    ecdf_contrasts.to_csv(output_directory / "conditional_ecdf_contrasts_bootstrap.csv", index=False)

    association_summary_frames: list[pd.DataFrame] = []
    association_draw_frames: list[pd.DataFrame] = []
    association_bootstrap_diagnostics: dict[str, Mapping[str, int]] = {}
    for target, aggregates in target_association_aggregates.items():
        summary, draws, diagnostics = _association_block_bootstrap(
            aggregates,
            replicates=int(arguments.bootstrap_replicates),
            seed=int(arguments.seed) + (105 if target == "lambda5" else 106),
        )
        association_summary_frames.append(summary)
        association_draw_frames.append(draws)
        association_bootstrap_diagnostics[target] = diagnostics
    association_summary = pd.concat(association_summary_frames, ignore_index=True)
    association_diagnostics = pd.concat(association_diagnostic_frames, ignore_index=True)
    association_draws = pd.concat(association_draw_frames, ignore_index=True)
    association_summary.to_csv(output_directory / "conditional_rank_associations.csv", index=False)
    association_diagnostics.to_csv(output_directory / "conditional_rank_association_diagnostics.csv", index=False)
    association_draws.to_csv(output_directory / "conditional_rank_association_bootstrap.csv", index=False)

    for target in ("lambda5", "lambda6"):
        _plot_ecdf(point_ecdf, ecdf_contrasts, target, figures_directory / f"conditional_ecdf_{target}.png")
    _plot_selected_pairs(pairs_frame, prepared, selected_pair_indices, figures_directory / "illustrative_conditional_pairs.png")
    _plot_tail_contrasts(contrasts, figures_directory / "conditional_tail_contrasts.png")
    _plot_conditional_rank_associations(association_summary, figures_directory / "conditional_rank_associations.png")

    _json_dump(output_directory / "source_hashes.json", source_hashes)
    _write_report(
        output_directory / "RUN_REPORT.md",
        source_hashes=source_hashes,
        support=support_frame,
        residual_summary=residual_frame,
        association_summary=association_summary,
        association_diagnostics=association_diagnostics,
        contrasts=contrasts,
        pair_rows=pairs_frame,
        bootstrap_diagnostics={
            "residual_group_contrasts": bootstrap_diagnostics,
            "fixed_rank_conditional_associations": association_bootstrap_diagnostics,
        },
    )
    _write_hash_manifest(output_directory)
    return {
        "status": "PASS",
        "output_directory": str(output_directory),
        "datasets": sorted(prepared),
        "bootstrap_replicates": int(arguments.bootstrap_replicates),
        "parity_rows": int(len(parity_frame)),
        "residual_rows": int(len(values_frame)),
        "conditional_association_rows": int(len(association_summary)),
        "conditional_association_bootstrap": association_bootstrap_diagnostics,
        "output_manifest": str(output_directory / "sha256_manifest.csv"),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Post-process Section 4.5 uniform windows without modifying a completed formal run."
    )
    parser.add_argument("run_dir", type=Path, help="Completed formal run containing direct_analysis_values.csv")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory outside run_dir for derivative evidence")
    parser.add_argument(
        "--common-config",
        type=Path,
        default=None,
        help="Optional byte-identical copy of the run's frozen common protocol; defaults to the run-local copy.",
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261012)
    parser.add_argument(
        "--block-shape",
        nargs=3,
        type=int,
        default=None,
        metavar=("X", "Y", "Z"),
        help="Must match the frozen formal primary block shape when supplied.",
    )
    parser.add_argument("--pair-caliper-iqr", type=float, default=1.0)
    parser.add_argument("--parity-sample-size", type=int, default=256)
    parser.add_argument("--parity-rtol", type=float, default=5.0e-6)
    parser.add_argument("--parity-atol", type=float, default=5.0e-8)
    arguments = parser.parse_args()
    if arguments.bootstrap_replicates <= 0 or arguments.parity_sample_size <= 0:
        parser.error("bootstrap and parity sample counts must be positive")
    if arguments.block_shape is not None and any(int(value) <= 0 for value in arguments.block_shape):
        parser.error("block dimensions must be positive")
    if arguments.pair_caliper_iqr <= 0.0:
        parser.error("pair caliper must be positive")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    try:
        result = run_analysis(arguments)
    except Exception as error:  # CLI boundary: preserve the failure reason in stderr.
        print(f"conditional_high_order_window_analysis: FAIL: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
