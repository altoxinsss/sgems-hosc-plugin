"""Numerical gates between the formal Python reference and ``movewindow_cli``.

The real-input gate retains the non-cubic 5 x 5 x 3 parity check.  A separate
fixed-fixture gate covers three window sizes, all formal weighting schemes,
orders one through six, and the frozen sensitivity settings.  Both archive
their CLI artifacts in self-contained JSON/CSV evidence bundles.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from cli_runner import (
    ATOL,
    RTOL,
    CUMULANT_CENTERING_VERSION,
    compare_maps,
    run_cli_case,
    sha256,
    validate_cli_window_alignment,
)
from grid_io import GridData, read_grid, write_map
from weighting import WeightConfig, compute_weighted_statistics, generate_weights
from window_support import extract_complete_windows


FORMULA_VERSION = "2026-10-10.v2"
FORMAL_WINDOW = (5, 5, 3)
FORMAL_ORDERS = (3, 4, 5, 6)
FORMAL_WEIGHT_CONFIGS: tuple[dict[str, Any], ...] = (
    {"mode": "uniform"},
    {"mode": "coordinate_gaussian", "bandwidth": 1.0},
    {
        "mode": "soft_high_porosity",
        "high_porosity_threshold": 0.25,
        "soft_transition_width": 0.02,
        "soft_weight_strength": 0.8,
    },
)
FORMAL_WEIGHT_MODES = frozenset({"uniform", "coordinate_gaussian", "soft_high_porosity"})
FIXTURE_WINDOWS = ((3, 3, 3), (5, 5, 3), (7, 7, 3))
FIXTURE_ORDERS = (1, 2, 3, 4, 5, 6)
FIXTURE_SENSITIVITY_CONFIGS: tuple[dict[str, Any], ...] = (
    {"id": "coordinate_gaussian_sigma_0.75", "mode": "coordinate_gaussian", "bandwidth": 0.75},
    {"id": "coordinate_gaussian_sigma_1.25", "mode": "coordinate_gaussian", "bandwidth": 1.25},
    {
        "id": "soft_high_porosity_width_0.01",
        "mode": "soft_high_porosity",
        "high_porosity_threshold": 0.25,
        "soft_transition_width": 0.01,
        "soft_weight_strength": 0.8,
    },
    {
        "id": "soft_high_porosity_width_0.04",
        "mode": "soft_high_porosity",
        "high_porosity_threshold": 0.25,
        "soft_transition_width": 0.04,
        "soft_weight_strength": 0.8,
    },
    {
        "id": "soft_high_porosity_strength_0.60",
        "mode": "soft_high_porosity",
        "high_porosity_threshold": 0.25,
        "soft_transition_width": 0.02,
        "soft_weight_strength": 0.6,
    },
    {
        "id": "soft_high_porosity_strength_1.00",
        "mode": "soft_high_porosity",
        "high_porosity_threshold": 0.25,
        "soft_transition_width": 0.02,
        "soft_weight_strength": 1.0,
    },
)


def _centering_metadata(
    metadata: Mapping[str, str],
    route: str,
    order: int,
) -> tuple[str | None, str | None, bool]:
    expected = "none" if int(order) == 1 else (
        "weighted_window_mean" if route == "weighted" else "window_mean"
    )
    centering = metadata.get("cumulant_centering")
    version = metadata.get("cumulant_centering_version")
    return centering, version, bool(
        centering == expected and version == CUMULANT_CENTERING_VERSION
    )


def _implementation_hashes() -> dict[str, str]:
    """Hash every source file that defines this numerical comparison."""
    experiment_dir = Path(__file__).resolve().parent
    repository_root = experiment_dir.parents[1]
    files = (
        Path(__file__),
        experiment_dir / "cli_runner.py",
        experiment_dir / "grid_io.py",
        experiment_dir / "statistics.py",
        experiment_dir / "weighting.py",
        experiment_dir / "window_support.py",
        repository_root / "src" / "window" / "main.cpp",
        repository_root / "vendor" / "weighted_moment_core" / "WeightedMomentCore.cpp",
        repository_root / "vendor" / "weighted_moment_core" / "WeightedMomentCore.h",
    )
    return {str(path.relative_to(repository_root)).replace("\\", "/"): sha256(path) for path in files}


@dataclass(frozen=True)
class ParityGateResult:
    """Locations and aggregate status of one immutable parity-gate execution."""

    passed: bool
    report_path: Path
    summary_path: Path
    max_abs_error: float
    case_count: int


@dataclass(frozen=True)
class FixtureParityGateResult:
    """Evidence locations for the exhaustive fixed-fixture C++ parity matrix."""

    passed: bool
    report_path: Path
    summary_path: Path
    max_abs_error: float
    case_count: int


@dataclass(frozen=True)
class UniformRouteAuditResult:
    """Evidence locations for the AVL/Bell/all-one-weighted uniform audit."""

    passed: bool
    report_path: Path
    summary_path: Path
    max_abs_error: float
    case_count: int


def default_movewindow_executable() -> Path:
    """Return the Release executable used by the checked-in test convention."""
    repository_root = Path(__file__).resolve().parents[2]
    default_build = repository_root / "build" / "Release" if os.name == "nt" else repository_root / "build"
    build_root = Path(os.environ.get("WITHOUT_SGEMS_BUILD_DIR", default_build))
    return build_root / ("movewindow_cli.exe" if os.name == "nt" else "movewindow_cli")


def _metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def _artifact_hashes(result) -> dict[str, str]:
    artifacts = {
        "ini_sha256": sha256(result.config_path),
        "output_sha256": sha256(result.output_path),
        "metadata_sha256": sha256(result.metadata_path),
        "stdout_sha256": sha256(result.stdout_path),
        "stderr_sha256": sha256(result.stderr_path),
    }
    if result.sidecar_path is not None:
        artifacts["window_support_sha256"] = sha256(result.sidecar_path)
    return artifacts


def _validate_protocol(
    grid: GridData,
    window_size: tuple[int, int, int],
    orders: Sequence[int],
    weight_configs: Sequence[Mapping[str, Any]],
) -> tuple[list[int], list[dict[str, Any]]]:
    if tuple(window_size) != FORMAL_WINDOW:
        raise ValueError(f"The formal C++/Python gate requires the non-cubic {FORMAL_WINDOW} window")
    if grid.values.shape[0] < window_size[2] or grid.values.shape[1] < window_size[1] or grid.values.shape[2] < window_size[0]:
        raise ValueError(f"Grid shape {grid.values.shape} does not contain a full {window_size} window")
    resolved_orders = [int(order) for order in orders]
    if set(resolved_orders) != set(FORMAL_ORDERS) or len(resolved_orders) != len(FORMAL_ORDERS):
        raise ValueError(f"The formal C++/Python gate requires orders {FORMAL_ORDERS}")
    resolved_weights = [dict(item) for item in weight_configs]
    modes = [WeightConfig.from_mapping(item).active_mode() for item in resolved_weights]
    if set(modes) != FORMAL_WEIGHT_MODES or len(modes) != len(FORMAL_WEIGHT_MODES):
        raise ValueError("The formal C++/Python gate requires one uniform, coordinate_gaussian, and soft_high_porosity case")
    return resolved_orders, resolved_weights


def _reference_maps(
    grid: GridData,
    window_size: tuple[int, int, int],
    move_step: tuple[int, int, int],
    weight_mapping: Mapping[str, Any],
    high_threshold: float,
    variance_epsilon: float,
    source_spacing: tuple[float, float, float],
    max_order: int = 6,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, float]]:
    batch = extract_complete_windows(grid, window_size, move_step)
    configuration = WeightConfig.from_mapping(weight_mapping)
    configuration.validate()
    weights = generate_weights(
        batch.values,
        configuration,
        window_size=window_size,
        source_spacing=source_spacing,
    )
    statistics = compute_weighted_statistics(
        batch,
        weights,
        high_threshold,
        variance_epsilon,
        cumulant_order=max_order,
    )
    diagnostics = {
        "window_count": int(len(batch.values)),
        "weight_sum_min": float(np.nanmin(statistics.weight_sum)),
        "weight_sum_max": float(np.nanmax(statistics.weight_sum)),
        "ess_min": float(np.nanmin(statistics.effective_sample_size)),
        "ess_max": float(np.nanmax(statistics.effective_sample_size)),
    }
    return batch.coordinates.copy(), statistics.statistics.as_dict(), diagnostics


def _reference_case(
    grid: GridData,
    window_size: tuple[int, int, int],
    move_step: tuple[int, int, int],
    weight_mapping: Mapping[str, Any],
    high_threshold: float,
    variance_epsilon: float,
    source_spacing: tuple[float, float, float],
    max_order: int,
):
    """Return the Python reference windows, weights, and statistics for one case."""

    batch = extract_complete_windows(grid, window_size, move_step)
    configuration = WeightConfig.from_mapping(weight_mapping)
    configuration.validate()
    weights = generate_weights(
        batch.values,
        configuration,
        window_size=window_size,
        source_spacing=source_spacing,
    )
    statistics = compute_weighted_statistics(
        batch,
        weights,
        high_threshold,
        variance_epsilon,
        cumulant_order=max_order,
    )
    return batch, weights, statistics


def _reference_cumulant(statistics, order: int) -> np.ndarray:
    fields = statistics.statistics.as_dict()
    if order == 1:
        return np.asarray(fields["mean"], dtype=float)
    if order == 2:
        return np.asarray(fields["cumulant2"], dtype=float)
    return np.asarray(fields[f"cumulant{order}"], dtype=float)


def _validate_sidecar_against_reference(result, batch, weights: np.ndarray) -> dict[str, float | int | str | bool]:
    """Check sidecar centers and weight diagnostics against the Python reference."""

    if result.sidecar_frame is None:
        raise RuntimeError("C++ parity case did not produce a window-support sidecar")
    alignment = validate_cli_window_alignment(
        result.frame,
        result.sidecar_frame,
        batch.coordinates,
        require_valid_support=True,
    )
    expected = batch.coordinates.copy()
    expected["weight_sum"] = np.sum(weights, axis=1)
    expected["weight_squared_sum"] = np.sum(weights**2, axis=1)
    expected["effective_sample_size"] = expected["weight_sum"] ** 2 / expected["weight_squared_sum"]
    expected["sample_count"] = int(weights.shape[1])
    observed = result.sidecar_frame.loc[result.sidecar_frame["map_included"]].merge(
        expected,
        on=["x", "y", "z"],
        how="inner",
        suffixes=("_cpp", "_python"),
    )
    if len(observed) != len(expected):
        raise RuntimeError("C++ sidecar diagnostics could not be aligned to the Python windows")
    errors: dict[str, float] = {}
    for name in ("weight_sum", "weight_squared_sum", "effective_sample_size"):
        difference = np.abs(observed[f"{name}_cpp"].to_numpy(float) - observed[f"{name}_python"].to_numpy(float))
        maximum = float(np.max(difference)) if len(difference) else 0.0
        if not np.allclose(observed[f"{name}_cpp"], observed[f"{name}_python"], rtol=RTOL, atol=ATOL):
            raise RuntimeError(f"C++ sidecar {name} disagrees with Python reference; max_abs_error={maximum}")
        errors[f"{name}_max_abs_error"] = maximum
    if not np.array_equal(observed["sample_count_cpp"].to_numpy(int), observed["sample_count_python"].to_numpy(int)):
        raise RuntimeError("C++ sidecar sample_count disagrees with the Python window support")
    return {
        "passed": True,
        "coordinate_validation": str(alignment["status"]),
        "sidecar_node_count": int(len(observed)),
        **errors,
    }


def run_cpp_python_parity_gate(
    executable: str | Path,
    grid: GridData,
    dataset_id: str,
    output_dir: str | Path,
    *,
    source_spacing: tuple[float, float, float] = (1.0, 1.0, 1.0),
    move_step: tuple[int, int, int] = (1, 1, 1),
    weight_configs: Sequence[Mapping[str, Any]] = FORMAL_WEIGHT_CONFIGS,
    orders: Sequence[int] = FORMAL_ORDERS,
    high_threshold: float = 0.25,
    variance_epsilon: float = 1.0e-12,
    timeout_seconds: float = 600.0,
    atol: float = ATOL,
    rtol: float = RTOL,
    fail_on_error: bool = True,
) -> ParityGateResult:
    """Run and archive the mandatory 3-weight x 4-order CLI parity matrix.

    The JSON report is written even when a case fails.  With ``fail_on_error``
    enabled, the caller receives an exception after those diagnostic artifacts
    exist, so a formal experiment cannot proceed with an unverified C++ core.
    """
    executable_path = Path(executable).resolve()
    if not executable_path.is_file():
        raise FileNotFoundError(f"movewindow_cli executable does not exist: {executable_path}")
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    resolved_orders, resolved_weights = _validate_protocol(grid, FORMAL_WINDOW, orders, weight_configs)
    source_spacing = tuple(float(value) for value in source_spacing)
    if len(source_spacing) != 3 or any(value <= 0.0 for value in source_spacing):
        raise ValueError("source_spacing must contain three positive values")
    if tuple(int(value) for value in move_step) != (1, 1, 1):
        raise ValueError("The formal parity gate requires unit move steps")

    cases: list[dict[str, Any]] = []
    for weight_mapping in resolved_weights:
        configuration = WeightConfig.from_mapping(weight_mapping)
        configuration.validate()
        mode = configuration.active_mode()
        coordinates, reference_statistics, diagnostics = _reference_maps(
            grid,
            FORMAL_WINDOW,
            move_step,
            weight_mapping,
            high_threshold,
            variance_epsilon,
            source_spacing,
        )
        for order in resolved_orders:
            case: dict[str, Any] = {
                "dataset_id": dataset_id,
                "weighting": mode,
                "order": int(order),
                "weight_config": dict(weight_mapping),
                "diagnostics": diagnostics,
                "passed": False,
            }
            try:
                expected = coordinates.copy()
                expected["value"] = np.asarray(reference_statistics[f"cumulant{order}"], dtype=float)
                result = run_cli_case(
                    executable_path,
                    grid,
                    f"{dataset_id}_{mode}",
                    FORMAL_WINDOW,
                    move_step,
                    "weighted",
                    dict(weight_mapping),
                    source_spacing,
                    output_path,
                    timeout_seconds,
                    order,
                )
                comparison = compare_maps(expected, result.frame, atol=atol, rtol=rtol)
                if result.sidecar_frame is None:
                    raise RuntimeError("C++ parity case did not produce a window-support sidecar")
                sidecar_alignment = validate_cli_window_alignment(
                    result.frame,
                    result.sidecar_frame,
                    coordinates,
                    require_valid_support=True,
                )
                metadata = _metadata(result.metadata_path)
                formula_matches = metadata.get("weighting_formula_version") == FORMULA_VERSION
                centering, centering_version, centering_matches = _centering_metadata(
                    metadata, "weighted", order
                )
                case.update(
                    comparison,
                    formula_version=metadata.get("weighting_formula_version"),
                    formula_version_matches=bool(formula_matches),
                    cumulant_centering=centering,
                    cumulant_centering_version=centering_version,
                    cumulant_centering_matches=centering_matches,
                    cli_config=str(result.config_path),
                    cli_output=str(result.output_path),
                    cli_metadata=str(result.metadata_path),
                    sidecar_coordinate_validation=sidecar_alignment["status"],
                    sidecar_path=str(result.sidecar_path) if result.sidecar_path is not None else None,
                    artifacts=_artifact_hashes(result),
                )
                case["passed"] = bool(
                    comparison["passed"] and formula_matches and centering_matches
                )
                if not formula_matches:
                    case["error"] = (
                        f"Expected weighting_formula_version={FORMULA_VERSION}; "
                        f"received {metadata.get('weighting_formula_version')!r}"
                    )
                elif not centering_matches:
                    case["error"] = (
                        f"Expected cumulant centering version={CUMULANT_CENTERING_VERSION!r} "
                        f"with weighted_window_mean; received version={centering_version!r}, "
                        f"mode={centering!r}"
                    )
            except Exception as error:  # Preserve the case record for the evidence bundle.
                case["error"] = f"{type(error).__name__}: {error}"
            cases.append(case)

    passed = bool(len(cases) == len(FORMAL_WEIGHT_MODES) * len(FORMAL_ORDERS) and all(case["passed"] for case in cases))
    finite_errors = [float(case["max_abs_error"]) for case in cases if np.isfinite(case.get("max_abs_error", np.nan))]
    max_abs_error = max(finite_errors) if finite_errors else None
    summary_columns = [
        "dataset_id", "weighting", "order", "passed", "count", "failed_count",
        "max_abs_error", "max_relative_error", "formula_version", "formula_version_matches",
        "cumulant_centering", "cumulant_centering_version", "cumulant_centering_matches", "error",
    ]
    summary = pd.DataFrame(cases)
    for column in summary_columns:
        if column not in summary:
            summary[column] = np.nan
    summary_path = output_path / "cpp_python_parity_gate.csv"
    summary.loc[:, summary_columns].to_csv(summary_path, index=False)
    report = {
        "gate": "cpp_python_parity_5x5x3_orders_3_to_6",
        "status": "PASS" if passed else "FAIL",
        "formula_version": FORMULA_VERSION,
        "window_size": list(FORMAL_WINDOW),
        "move_step": list(move_step),
        "orders": resolved_orders,
        "input_grid": str(grid.source_path),
        "input_grid_sha256": sha256(grid.source_path),
        "executable": str(executable_path),
        "executable_sha256": sha256(executable_path),
        "implementation_sha256": _implementation_hashes(),
        "source_spacing": list(source_spacing),
        "absolute_tolerance": float(atol),
        "relative_tolerance": float(rtol),
        "case_count": len(cases),
        "max_abs_error": max_abs_error,
        "cases": cases,
    }
    report_path = output_path / "cpp_python_parity_gate.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    gate_result = ParityGateResult(
        passed,
        report_path,
        summary_path,
        float(max_abs_error) if max_abs_error is not None else float("nan"),
        len(cases),
    )
    if not passed and fail_on_error:
        raise RuntimeError(f"C++/Python parity gate failed; inspect {report_path}")
    return gate_result


def _write_fixture_grid(path: Path) -> GridData:
    """Create a deterministic non-Gaussian 9 x 9 x 5 parity fixture."""

    coordinates = pd.DataFrame(
        [(x, y, z) for z in range(5) for y in range(9) for x in range(9)],
        columns=["x", "y", "z"],
    )
    x = coordinates["x"].to_numpy(float)
    y = coordinates["y"].to_numpy(float)
    z = coordinates["z"].to_numpy(float)
    values = (
        0.14
        + 0.010 * x
        + 0.008 * y
        + 0.016 * z
        + 0.0014 * x * y
        - 0.0009 * y * z
        + 0.0005 * x**2
        + 0.006 * np.sin(0.7 * x + 0.4 * y - 0.3 * z)
    )
    write_map(path, coordinates, values, "value", title="Fixed C++/Python parity fixture")
    return read_grid(path, "value")


def _run_fixture_case(
    executable: Path,
    grid: GridData,
    output_dir: Path,
    *,
    window_size: tuple[int, int, int],
    weight_mapping: Mapping[str, Any],
    parameter_role: str,
    order: int,
    source_spacing: tuple[float, float, float],
    high_threshold: float,
    variance_epsilon: float,
    timeout_seconds: float,
    atol: float,
    rtol: float,
) -> dict[str, Any]:
    configuration = WeightConfig.from_mapping(weight_mapping)
    configuration.validate()
    mode = configuration.active_mode()
    if mode not in FORMAL_WEIGHT_MODES:
        raise ValueError(f"Fixture parity only supports formal weights, got {mode}")
    parameter_id = str(weight_mapping.get("id", mode))
    route = "avl" if mode == "uniform" else "weighted"
    batch, weights, statistics = _reference_case(
        grid,
        window_size,
        (1, 1, 1),
        weight_mapping,
        high_threshold,
        variance_epsilon,
        source_spacing,
        max(FIXTURE_ORDERS),
    )
    expected = batch.coordinates.copy()
    expected["value"] = _reference_cumulant(statistics, order)
    case_id = (
        f"fixture_{parameter_role}_{parameter_id}_"
        f"w{'x'.join(map(str, window_size))}_kappa{order}_{route}"
    )
    result = run_cli_case(
        executable,
        grid,
        "fixed_fixture",
        window_size,
        (1, 1, 1),
        route,
        dict(weight_mapping),
        source_spacing,
        output_dir,
        timeout_seconds,
        order,
        case_id=case_id,
        target_property=f"kappa{order}",
        require_window_support_sidecar=True,
    )
    comparison = compare_maps(expected, result.frame, atol=atol, rtol=rtol)
    sidecar = _validate_sidecar_against_reference(result, batch, weights)
    metadata = _metadata(result.metadata_path)
    formula_matches = metadata.get("weighting_formula_version") == FORMULA_VERSION
    centering, centering_version, centering_matches = _centering_metadata(
        metadata, route, order
    )
    return {
        "parameter_role": parameter_role,
        "parameter_id": parameter_id,
        "weighting": mode,
        "weight_config": dict(weight_mapping),
        "window_size": list(window_size),
        "order": int(order),
        "algorithm": "AVL" if route == "avl" else "Weighted moment",
        **comparison,
        "formula_version": metadata.get("weighting_formula_version"),
        "formula_version_matches": bool(formula_matches),
        "cumulant_centering": centering,
        "cumulant_centering_version": centering_version,
        "cumulant_centering_matches": centering_matches,
        "sidecar": sidecar,
        "cli_config": str(result.config_path),
        "cli_output": str(result.output_path),
        "cli_metadata": str(result.metadata_path),
        "cli_sidecar": str(result.sidecar_path) if result.sidecar_path is not None else None,
        "artifacts": _artifact_hashes(result),
        "passed": bool(
            comparison["passed"] and formula_matches and centering_matches and sidecar["passed"]
        ),
    }


def run_cpp_python_fixture_matrix_gate(
    executable: str | Path,
    output_dir: str | Path,
    *,
    source_spacing: tuple[float, float, float] = (2.0, 3.0, 4.0),
    high_threshold: float = 0.25,
    variance_epsilon: float = 1.0e-12,
    timeout_seconds: float = 120.0,
    atol: float = ATOL,
    rtol: float = RTOL,
    fail_on_error: bool = True,
) -> FixtureParityGateResult:
    """Run the complete small-grid C++/Python parity matrix.

    This gate covers the three formal windows, all three primary weights, and
    orders one through six (54 cases).  It also runs one sixth-order fixed
    fixture for each of the six frozen sensitivity parameter settings.  It is
    deliberately separate from the real-data 5x5x3 gate so formal runners can
    archive both a broad deterministic regression gate and a real-input gate.
    """

    executable_path = Path(executable).resolve()
    if not executable_path.is_file():
        raise FileNotFoundError(f"movewindow_cli executable does not exist: {executable_path}")
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    spacing = tuple(float(value) for value in source_spacing)
    if len(spacing) != 3 or any(value <= 0.0 for value in spacing):
        raise ValueError("source_spacing must contain three positive values")
    grid = _write_fixture_grid(output_path / "cpp_python_fixture_9x9x5.grid")
    cases: list[dict[str, Any]] = []

    for window_size in FIXTURE_WINDOWS:
        for weight_mapping in FORMAL_WEIGHT_CONFIGS:
            for order in FIXTURE_ORDERS:
                case: dict[str, Any] = {
                    "parameter_role": "primary",
                    "weighting": WeightConfig.from_mapping(weight_mapping).active_mode(),
                    "window_size": list(window_size),
                    "order": int(order),
                    "passed": False,
                }
                try:
                    case.update(
                        _run_fixture_case(
                            executable_path,
                            grid,
                            output_path,
                            window_size=window_size,
                            weight_mapping=weight_mapping,
                            parameter_role="primary",
                            order=order,
                            source_spacing=spacing,
                            high_threshold=high_threshold,
                            variance_epsilon=variance_epsilon,
                            timeout_seconds=timeout_seconds,
                            atol=atol,
                            rtol=rtol,
                        )
                    )
                except Exception as error:
                    case["error"] = f"{type(error).__name__}: {error}"
                cases.append(case)

    for weight_mapping in FIXTURE_SENSITIVITY_CONFIGS:
        case = {
            "parameter_role": "sensitivity",
            "weighting": WeightConfig.from_mapping(weight_mapping).active_mode(),
            "window_size": [5, 5, 3],
            "order": 6,
            "passed": False,
        }
        try:
            case.update(
                _run_fixture_case(
                    executable_path,
                    grid,
                    output_path,
                    window_size=(5, 5, 3),
                    weight_mapping=weight_mapping,
                    parameter_role="sensitivity",
                    order=6,
                    source_spacing=spacing,
                    high_threshold=high_threshold,
                    variance_epsilon=variance_epsilon,
                    timeout_seconds=timeout_seconds,
                    atol=atol,
                    rtol=rtol,
                )
            )
        except Exception as error:
            case["error"] = f"{type(error).__name__}: {error}"
        cases.append(case)

    expected_case_count = len(FIXTURE_WINDOWS) * len(FORMAL_WEIGHT_CONFIGS) * len(FIXTURE_ORDERS) + len(FIXTURE_SENSITIVITY_CONFIGS)
    passed = bool(len(cases) == expected_case_count and all(case["passed"] for case in cases))
    finite_errors = [float(case["max_abs_error"]) for case in cases if np.isfinite(case.get("max_abs_error", np.nan))]
    max_abs_error = max(finite_errors) if finite_errors else None
    summary_columns = [
        "parameter_role", "parameter_id", "weighting", "window_size", "order", "algorithm", "passed",
        "count", "failed_count", "max_abs_error", "max_relative_error", "formula_version",
        "formula_version_matches", "error",
    ]
    summary = pd.DataFrame(cases)
    for column in summary_columns:
        if column not in summary:
            summary[column] = np.nan
    summary_path = output_path / "cpp_python_fixture_matrix_gate.csv"
    summary.loc[:, summary_columns].to_csv(summary_path, index=False)
    report = {
        "gate": "cpp_python_fixture_matrix_windows_3_5_7_orders_1_to_6",
        "status": "PASS" if passed else "FAIL",
        "formula_version": FORMULA_VERSION,
        "fixture_grid": str(grid.source_path),
        "fixture_grid_sha256": sha256(grid.source_path),
        "fixture_shape_zyx": list(grid.values.shape),
        "executable": str(executable_path),
        "executable_sha256": sha256(executable_path),
        "implementation_sha256": _implementation_hashes(),
        "primary_windows": [list(window) for window in FIXTURE_WINDOWS],
        "primary_orders": list(FIXTURE_ORDERS),
        "primary_weight_configs": list(FORMAL_WEIGHT_CONFIGS),
        "sensitivity_weight_configs": list(FIXTURE_SENSITIVITY_CONFIGS),
        "source_spacing": list(spacing),
        "absolute_tolerance": float(atol),
        "relative_tolerance": float(rtol),
        "expected_case_count": expected_case_count,
        "case_count": len(cases),
        "max_abs_error": max_abs_error,
        "cases": cases,
    }
    report_path = output_path / "cpp_python_fixture_matrix_gate.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    result = FixtureParityGateResult(
        passed,
        report_path,
        summary_path,
        float(max_abs_error) if max_abs_error is not None else float("nan"),
        len(cases),
    )
    if not passed and fail_on_error:
        raise RuntimeError(f"C++/Python fixture matrix parity gate failed; inspect {report_path}")
    return result


def run_uniform_three_route_audit(
    executable: str | Path,
    output_dir: str | Path,
    *,
    source_spacing: tuple[float, float, float] = (2.0, 3.0, 4.0),
    timeout_seconds: float = 120.0,
    atol: float = ATOL,
    rtol: float = RTOL,
    fail_on_error: bool = True,
) -> UniformRouteAuditResult:
    """Compare the three C++ routes on one frozen non-Gaussian uniform case.

    Uniform production maps use AVL.  This audit additionally evaluates the
    same fixed windows through Bell and the weighted-moment route with exactly
    one weight at every sample.  It verifies the values and the all-one weight
    diagnostics for every cumulant order one through six.
    """
    executable_path = Path(executable).resolve()
    if not executable_path.is_file():
        raise FileNotFoundError(f"movewindow_cli executable does not exist: {executable_path}")
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    spacing = tuple(float(value) for value in source_spacing)
    if len(spacing) != 3 or any(value <= 0.0 for value in spacing):
        raise ValueError("source_spacing must contain three positive values")

    grid = _write_fixture_grid(output_path / "uniform_three_route_fixture_9x9x5.grid")
    uniform = {"mode": "uniform"}
    window = (5, 5, 3)
    cases: list[dict[str, Any]] = []
    for order in FIXTURE_ORDERS:
        case: dict[str, Any] = {"order": int(order), "passed": False}
        try:
            routes = {
                route: run_cli_case(
                    executable_path,
                    grid,
                    "uniform_three_route_fixture",
                    window,
                    (1, 1, 1),
                    route,
                    uniform,
                    spacing,
                    output_path,
                    timeout_seconds,
                    order,
                    case_id=f"uniform_three_route_kappa{order}_{route}",
                    target_property=f"kappa{order}_{route}",
                    require_window_support_sidecar=True,
                )
                for route in ("avl", "bell", "weighted")
            }
            comparisons = {
                "bell_minus_avl": compare_maps(routes["avl"].frame, routes["bell"].frame, atol=atol, rtol=rtol),
                "weighted_minus_avl": compare_maps(routes["avl"].frame, routes["weighted"].frame, atol=atol, rtol=rtol),
            }
            weighted_sidecar = routes["weighted"].sidecar_frame
            if weighted_sidecar is None:
                raise RuntimeError("All-one weighted route did not produce a support sidecar")
            included = weighted_sidecar.loc[weighted_sidecar["map_included"].astype(bool)].copy()
            expected_count = included["sample_count"].to_numpy(float)
            all_one = bool(
                np.allclose(included["weight_sum"].to_numpy(float), expected_count, atol=atol, rtol=rtol)
                and np.allclose(included["weight_squared_sum"].to_numpy(float), expected_count, atol=atol, rtol=rtol)
                and np.allclose(included["effective_sample_size"].to_numpy(float), expected_count, atol=atol, rtol=rtol)
            )
            route_metadata = {
                route: _metadata(result.metadata_path) for route, result in routes.items()
            }
            weighted_metadata = route_metadata["weighted"]
            formula_matches = weighted_metadata.get("weighting_formula_version") == FORMULA_VERSION
            centering_checks = {
                route: _centering_metadata(metadata, route, order)
                for route, metadata in route_metadata.items()
            }
            centering_matches = all(item[2] for item in centering_checks.values())
            max_abs_error = max(
                float(comparisons["bell_minus_avl"]["max_abs_error"]),
                float(comparisons["weighted_minus_avl"]["max_abs_error"]),
            )
            case.update(
                comparisons,
                all_one_weight_diagnostics=all_one,
                formula_version=weighted_metadata.get("weighting_formula_version"),
                formula_version_matches=formula_matches,
                cumulant_centering={route: item[0] for route, item in centering_checks.items()},
                cumulant_centering_version={route: item[1] for route, item in centering_checks.items()},
                cumulant_centering_matches=centering_matches,
                max_abs_error=max_abs_error,
                artifacts={route: _artifact_hashes(result) for route, result in routes.items()},
            )
            case["passed"] = bool(
                comparisons["bell_minus_avl"]["passed"]
                and comparisons["weighted_minus_avl"]["passed"]
                and all_one
                and formula_matches
                and centering_matches
            )
        except Exception as error:
            case["error"] = f"{type(error).__name__}: {error}"
        cases.append(case)

    passed = bool(len(cases) == len(FIXTURE_ORDERS) and all(case["passed"] for case in cases))
    finite_errors = [float(case["max_abs_error"]) for case in cases if np.isfinite(case.get("max_abs_error", np.nan))]
    max_abs_error = max(finite_errors) if finite_errors else None
    summary_columns = [
        "order", "passed", "max_abs_error", "all_one_weight_diagnostics",
        "formula_version", "formula_version_matches", "error",
    ]
    summary = pd.DataFrame(cases)
    for column in summary_columns:
        if column not in summary:
            summary[column] = np.nan
    summary_path = output_path / "uniform_three_route_audit.csv"
    summary.loc[:, summary_columns].to_csv(summary_path, index=False)
    report = {
        "gate": "uniform_avl_bell_all_one_weighted_orders_1_to_6",
        "status": "PASS" if passed else "FAIL",
        "fixture_grid": str(grid.source_path),
        "fixture_grid_sha256": sha256(grid.source_path),
        "window_size": list(window),
        "orders": list(FIXTURE_ORDERS),
        "source_spacing": list(spacing),
        "formula_version": FORMULA_VERSION,
        "executable": str(executable_path),
        "executable_sha256": sha256(executable_path),
        "absolute_tolerance": float(atol),
        "relative_tolerance": float(rtol),
        "case_count": len(cases),
        "max_abs_error": max_abs_error,
        "cases": cases,
    }
    report_path = output_path / "uniform_three_route_audit.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    result = UniformRouteAuditResult(
        passed,
        report_path,
        summary_path,
        float(max_abs_error) if max_abs_error is not None else float("nan"),
        len(cases),
    )
    if not passed and fail_on_error:
        raise RuntimeError(f"Uniform three-route audit failed; inspect {report_path}")
    return result


def _triple(raw: str, *, kind: type[int] | type[float]) -> tuple:
    values = tuple(kind(value.strip()) for value in raw.split(","))
    if len(values) != 3:
        raise argparse.ArgumentTypeError("Expected exactly three comma-separated values")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the formal C++/Python 5x5x3 parity gate")
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--property", default=None)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=default_movewindow_executable())
    parser.add_argument("--source-spacing", default="1,1,1")
    parser.add_argument("--timeout-seconds", type=float, default=600.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    spacing = _triple(args.source_spacing, kind=float)
    result = run_cpp_python_parity_gate(
        args.executable,
        read_grid(args.grid, args.property),
        args.dataset_id,
        args.output_dir,
        source_spacing=spacing,
        timeout_seconds=args.timeout_seconds,
    )
    print(result.report_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
