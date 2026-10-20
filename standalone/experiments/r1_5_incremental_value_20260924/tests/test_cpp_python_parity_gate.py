from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from cpp_python_parity_gate import (
    FIXTURE_ORDERS,
    FIXTURE_SENSITIVITY_CONFIGS,
    FIXTURE_WINDOWS,
    FORMAL_ORDERS,
    FORMAL_WEIGHT_CONFIGS,
    FORMULA_VERSION,
    run_cpp_python_fixture_matrix_gate,
    run_cpp_python_parity_gate,
    run_uniform_three_route_audit,
)
from grid_io import read_grid, write_map


def weighted_executable(repo_root: Path) -> Path:
    default_build = repo_root / "build" / "Release" if os.name == "nt" else repo_root / "build"
    return Path(os.environ.get("WITHOUT_SGEMS_BUILD_DIR", default_build)) / ("movewindow_cli.exe" if os.name == "nt" else "movewindow_cli")


def test_formal_cpp_python_parity_gate_covers_all_weights_and_orders(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    coordinates = pd.DataFrame(
        [(x, y, z) for z in range(4) for y in range(6) for x in range(7)],
        columns=["x", "y", "z"],
    )
    x, y, z = (coordinates[column].to_numpy(float) for column in ("x", "y", "z"))
    values = (
        0.12 + 0.011 * x + 0.017 * y + 0.023 * z
        + 0.0011 * x * y - 0.0007 * y * z + 0.0003 * x**2
    )
    source = tmp_path / "non_cubic_volume.grid"
    write_map(source, coordinates, values, "value")

    result = run_cpp_python_parity_gate(
        executable,
        read_grid(source, "value"),
        "synthetic_non_cubic",
        tmp_path / "parity_gate",
        source_spacing=(2.0, 3.0, 4.0),
        timeout_seconds=60.0,
    )

    assert result.passed
    assert result.case_count == 12
    assert result.max_abs_error <= 1.0e-12
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["window_size"] == [5, 5, 3]
    assert report["orders"] == list(FORMAL_ORDERS)
    assert len(report["cases"]) == 12
    assert {case["weighting"] for case in report["cases"]} == {
        "uniform", "coordinate_gaussian", "soft_high_porosity",
    }
    assert {case["formula_version"] for case in report["cases"]} == {FORMULA_VERSION}
    assert all(case["passed"] for case in report["cases"])
    assert result.summary_path.is_file()


def test_fixture_cpp_python_parity_gate_covers_orders_windows_primary_and_sensitivity_cases(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    result = run_cpp_python_fixture_matrix_gate(
        weighted_executable(repo_root),
        tmp_path / "fixture_matrix_gate",
        timeout_seconds=60.0,
    )

    expected_cases = len(FIXTURE_WINDOWS) * len(FORMAL_WEIGHT_CONFIGS) * len(FIXTURE_ORDERS) + len(FIXTURE_SENSITIVITY_CONFIGS)
    assert result.passed
    assert result.case_count == expected_cases
    assert result.max_abs_error <= 1.0e-12
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["primary_orders"] == list(FIXTURE_ORDERS)
    assert report["primary_windows"] == [list(window) for window in FIXTURE_WINDOWS]
    assert len(report["cases"]) == expected_cases
    assert all(case["passed"] for case in report["cases"])
    primary = [case for case in report["cases"] if case["parameter_role"] == "primary"]
    sensitivity = [case for case in report["cases"] if case["parameter_role"] == "sensitivity"]
    assert len(primary) == 54
    assert len(sensitivity) == len(FIXTURE_SENSITIVITY_CONFIGS)
    assert {case["order"] for case in primary} == set(FIXTURE_ORDERS)
    assert {tuple(case["window_size"]) for case in primary} == set(FIXTURE_WINDOWS)
    assert {case["parameter_id"] for case in sensitivity} == {item["id"] for item in FIXTURE_SENSITIVITY_CONFIGS}
    assert all(case["sidecar"]["coordinate_validation"] == "PASS" for case in report["cases"])


def test_uniform_three_route_audit_matches_avl_bell_and_all_one_weighted(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    result = run_uniform_three_route_audit(
        weighted_executable(repo_root),
        tmp_path / "uniform_three_route",
        timeout_seconds=60.0,
    )

    assert result.passed
    assert result.case_count == len(FIXTURE_ORDERS)
    assert result.max_abs_error <= 1.0e-12
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert all(case["passed"] for case in report["cases"])
    assert all(case["all_one_weight_diagnostics"] for case in report["cases"])
