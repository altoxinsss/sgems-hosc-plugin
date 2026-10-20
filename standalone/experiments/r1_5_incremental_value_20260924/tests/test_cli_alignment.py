import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cli_runner import (
    CppFeatureMapSpec,
    WINDOW_SUPPORT_COLUMNS,
    _experiment_weight_switches,
    build_cpp_feature_map_specs,
    compare_maps,
    run_cli_case,
    run_cpp_feature_matrix,
    validate_cpp_feature_map_protocol,
    validate_cli_window_alignment,
)
from grid_io import read_grid, write_map
from weighting import WeightConfig, compute_weighted_statistics, generate_weights
from window_support import extract_complete_windows


def weighted_executable(repo_root: Path) -> Path:
    default_build = repo_root / "build" / "Release" if os.name == "nt" else repo_root / "build"
    return Path(os.environ.get("WITHOUT_SGEMS_BUILD_DIR", default_build)) / ("movewindow_cli.exe" if os.name == "nt" else "movewindow_cli")


def metadata_value(path: Path, key: str) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(key + "="):
            return float(line.split("=", 1)[1])
    raise AssertionError(f"Missing metadata key {key} in {path}")


def test_bell_avl_alignment_on_fixture(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    fixture = repo_root / "tests" / "fixtures" / "window_3x3.grid"
    grid = read_grid(fixture, "value")
    weight = {"mode": "uniform", "min_weight": 0.0, "normalize": True}
    bell = run_cli_case(executable, grid, "fixture", (3, 3, 1), (1, 1, 0), "bell", weight, (1, 1, 1), tmp_path, 30)
    avl = run_cli_case(executable, grid, "fixture", (3, 3, 1), (1, 1, 0), "avl", weight, (1, 1, 1), tmp_path, 30)
    assert compare_maps(bell.frame, avl.frame)["passed"]
    assert avl.sidecar_path is not None and avl.sidecar_path.is_file()
    assert avl.sidecar_frame is not None
    assert tuple(avl.sidecar_frame.columns) == WINDOW_SUPPORT_COLUMNS
    assert (avl.sidecar_frame["map_included"] & avl.sidecar_frame["valid"]).all()
    assert np.allclose(avl.sidecar_frame["weight_sum"], 9.0)
    assert np.allclose(avl.sidecar_frame["weight_squared_sum"], 9.0)
    assert np.allclose(avl.sidecar_frame["effective_sample_size"], 9.0)


@pytest.mark.parametrize(
    ("weight", "expected_scheme"),
    [
        ({"mode": "uniform", "min_weight": 0.0, "normalize": True}, "uniform.v1"),
        (
            {"mode": "soft_high_porosity", "enable_soft_high_porosity": True, "high_porosity_threshold": 0.25,
             "soft_transition_width": 0.02, "soft_weight_strength": 0.8, "min_weight": 0.0, "normalize": True},
            "soft_high_porosity.v1",
        ),
    ],
)
def test_cpp_python_weighting_scheme_alignment(tmp_path: Path, weight: dict, expected_scheme: str):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    grid = read_grid(repo_root / "tests" / "fixtures" / "window_3x3.grid", "value")
    result = run_cli_case(executable, grid, "fixture", (3, 3, 1), (1, 1, 0), "weighted", weight, (1, 1, 1), tmp_path, 30)
    batch = extract_complete_windows(grid, (3, 3, 1), (1, 1, 0))
    weights = generate_weights(batch.values, WeightConfig.from_mapping(weight), window_size=(3, 3, 1))
    expected = batch.coordinates.copy()
    expected["value"] = compute_weighted_statistics(batch, weights, 0.25, 1.0e-15).statistics.cumulant4
    assert compare_maps(expected, result.frame)["passed"]
    metadata = result.metadata_path.read_text(encoding="utf-8")
    assert f"weighting_scheme_id={expected_scheme}" in metadata
    assert "weighting_formula_version=2026-10-10.v2" in metadata


def test_coordinate_gaussian_cpp_python_alignment_and_diagnostics_on_3d_volume(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    coordinates = pd.DataFrame([(x, y, z) for z in range(3) for y in range(3) for x in range(3)], columns=["x", "y", "z"])
    source = tmp_path / "volume_3x3x3.grid"
    values = 0.12 + 0.01 * coordinates["x"] + 0.02 * coordinates["y"] + 0.03 * coordinates["z"] + 0.001 * coordinates["x"] * coordinates["z"]
    write_map(source, coordinates, values, "value")
    grid = read_grid(source, "value")
    weight = {"mode": "coordinate_gaussian", "enable_coordinate_gaussian": True, "bandwidth": 2.0, "min_weight": 0.0, "normalize": True}
    spacing = (2.0, 3.0, 4.0)
    result = run_cli_case(executable, grid, "volume", (3, 3, 3), (1, 1, 1), "weighted", weight, spacing, tmp_path, 30)
    ini = result.config_path.read_text(encoding="utf-8")
    assert "coordinate_dimensions=3" in ini
    assert "weight_center_mode=window_center" in ini
    batch = extract_complete_windows(grid, (3, 3, 3), (1, 1, 1))
    weights = generate_weights(batch.values, WeightConfig.from_mapping(weight), window_size=(3, 3, 3), source_spacing=spacing)
    expected = batch.coordinates.copy()
    expected["value"] = compute_weighted_statistics(batch, weights, 0.25, 1.0e-15).statistics.cumulant4
    assert compare_maps(expected, result.frame)["passed"]
    weight_sum = float(np.sum(weights[0]))
    square_sum = float(np.sum(weights[0] ** 2))
    ess = weight_sum**2 / square_sum
    assert np.isclose(metadata_value(result.metadata_path, "weight_sum_mean"), weight_sum)
    assert np.isclose(metadata_value(result.metadata_path, "weight_square_sum_mean"), square_sum)
    assert np.isclose(metadata_value(result.metadata_path, "effective_sample_size_mean"), ess)
    assert result.sidecar_frame is not None
    support = result.sidecar_frame
    alignment = validate_cli_window_alignment(result.frame, support, batch.coordinates, require_valid_support=True)
    assert alignment["status"] == "PASS"
    assert np.allclose(support["weight_sum"].to_numpy(float), np.sum(weights, axis=1))
    assert np.allclose(support["weight_squared_sum"].to_numpy(float), np.sum(weights**2, axis=1))
    assert np.allclose(support["effective_sample_size"].to_numpy(float), np.sum(weights, axis=1) ** 2 / np.sum(weights**2, axis=1))


@pytest.mark.parametrize("order", [2, 3, 4, 5, 6])
@pytest.mark.parametrize(
    ("route", "weight", "expected_centering"),
    [
        ("avl", {"mode": "uniform"}, "window_mean"),
        ("bell", {"mode": "uniform"}, "window_mean"),
        ("weighted", {"mode": "uniform"}, "weighted_window_mean"),
        (
            "weighted",
            {"mode": "coordinate_gaussian", "bandwidth": 1.0},
            "weighted_window_mean",
        ),
        (
            "weighted",
            {
                "mode": "soft_high_porosity",
                "high_porosity_threshold": 0.25,
                "soft_transition_width": 0.02,
                "soft_weight_strength": 0.8,
            },
            "weighted_window_mean",
        ),
    ],
)
def test_high_order_cumulants_are_stable_for_near_constant_windows(
    tmp_path: Path,
    order: int,
    route: str,
    weight: dict,
    expected_centering: str,
) -> None:
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    coordinates = pd.DataFrame(
        [(x, y, z) for z in range(3) for y in range(3) for x in range(3)],
        columns=["x", "y", "z"],
    )
    # This is the z15_26 window that exposed raw-moment cancellation after
    # standardization during the formal run.
    values = np.asarray(
        [
            0.29545, 0.2945625, 0.2948875, 0.294225, 0.296225, 0.2953625,
            0.2944125, 0.295675, 0.2930625, 0.2988, 0.2967875, 0.29855,
            0.296025, 0.2970125, 0.29805, 0.2967, 0.2988, 0.297925,
            0.298575, 0.2986125, 0.2988, 0.29815, 0.2980375, 0.2965125,
            0.2975125, 0.2966375, 0.2971625,
        ],
        dtype=float,
    )
    source = tmp_path / f"near_constant_{route}_{order}.grid"
    write_map(source, coordinates, values, "value")
    grid = read_grid(source, "value")
    batch = extract_complete_windows(grid, (3, 3, 3), (1, 1, 1))
    resolved_weights = generate_weights(
        batch.values,
        WeightConfig.from_mapping(weight),
        window_size=(3, 3, 3),
        source_spacing=(50.0, 50.0, 2.0),
    )
    statistics = compute_weighted_statistics(
        batch, resolved_weights, 0.25, 1.0e-12, cumulant_order=6
    ).statistics.as_dict()
    expected = float(statistics[f"cumulant{order}"][0])
    variance = float(statistics["variance"][0])

    result = run_cli_case(
        executable,
        grid,
        "near_constant",
        (3, 3, 3),
        (1, 1, 1),
        route,
        weight,
        (50.0, 50.0, 2.0),
        tmp_path / f"run_{route}_{order}",
        30,
        order=order,
    )
    observed = float(result.frame["value"].iloc[0])
    standardized_error = abs(observed - expected) / variance ** (order / 2.0)
    assert standardized_error <= 1.0e-8
    ini = result.config_path.read_text(encoding="utf-8")
    assert "center_cumulant_input_for_order_ge_2=true" in ini
    metadata = result.metadata_path.read_text(encoding="utf-8")
    assert "center_cumulant_input_for_order_ge_2=true" in metadata
    assert f"cumulant_centering={expected_centering}" in metadata
    assert "cumulant_centering_version=order_ge_2.v1" in metadata


def test_cpp_feature_matrix_returns_authoritative_coordinate_validated_maps(tmp_path: Path):
    experiment_dir = Path(__file__).resolve().parents[1]
    repo_root = experiment_dir.parents[1]
    executable = weighted_executable(repo_root)
    coordinates = pd.DataFrame([(x, y, z) for z in range(3) for y in range(5) for x in range(5)], columns=["x", "y", "z"])
    source = tmp_path / "matrix_fixture.grid"
    values = 0.1 + 0.01 * coordinates["x"] + 0.02 * coordinates["y"] + 0.03 * coordinates["z"] + 0.002 * coordinates["x"] * coordinates["y"]
    write_map(source, coordinates, values, "value")
    grid = read_grid(source, "value")
    result = run_cpp_feature_matrix(
        executable,
        (
            CppFeatureMapSpec("fixture", grid, (3, 3, 3), 1, {"mode": "uniform"}),
            CppFeatureMapSpec("fixture", grid, (3, 3, 3), 1, {"mode": "coordinate_gaussian", "bandwidth": 1.0}),
        ),
        tmp_path / "matrix",
        timeout_seconds=30.0,
    )
    assert len(result.artifacts) == 2
    assert result.manifest_path.is_file()
    assert result.uniform_manifest_path.is_file()
    assert result.weighted_manifest_path.is_file()
    assert set(result.manifest["algorithm"]) == {"AVL", "Weighted moment"}
    assert set(result.manifest["cumulant_centering"]) == {"none"}
    assert set(result.manifest["cumulant_centering_version"]) == {"order_ge_2.v1"}
    assert result.manifest["authoritative_cplusplus"].all()
    assert set(result.manifest["coordinate_validation"]) == {"PASS"}
    for artifact in result.artifacts:
        assert artifact.authoritative_cplusplus
        assert artifact.cli_result.authoritative_cplusplus
        assert artifact.cli_result.sidecar_path is not None
        assert artifact.map_frame.columns.tolist() == ["x", "y", "z", "value"]
        assert not artifact.support_frame.empty
    uniform = next(item for item in result.artifacts if item.algorithm == "AVL")
    expected = extract_complete_windows(grid, (3, 3, 3), (1, 1, 1)).coordinates.copy()
    expected["value"] = np.mean(extract_complete_windows(grid, (3, 3, 3), (1, 1, 1)).values, axis=1)
    assert compare_maps(expected, uniform.map_frame)["passed"]


def test_cpp_feature_map_spec_builder_constructs_but_does_not_run_the_486_map_protocol(tmp_path: Path):
    coordinates = pd.DataFrame([(x, y, z) for z in range(3) for y in range(7) for x in range(7)], columns=["x", "y", "z"])
    source = tmp_path / "protocol_fixture.grid"
    write_map(source, coordinates, 0.15 + 0.01 * coordinates["x"] + 0.02 * coordinates["y"], "value")
    grid = read_grid(source, "value")
    primary = (
        {"id": "uniform", "mode": "uniform"},
        {"id": "coordinate_gaussian", "mode": "coordinate_gaussian", "bandwidth": 1.0},
        {"id": "soft_high_porosity", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 0.8},
    )
    sensitivity = (
        {"id": "coordinate_gaussian_sigma_0.75", "mode": "coordinate_gaussian", "bandwidth": 0.75},
        {"id": "coordinate_gaussian_sigma_1.25", "mode": "coordinate_gaussian", "bandwidth": 1.25},
        {"id": "soft_width_0.01", "mode": "soft_high_porosity", "soft_transition_width": 0.01, "soft_weight_strength": 0.8},
        {"id": "soft_width_0.04", "mode": "soft_high_porosity", "soft_transition_width": 0.04, "soft_weight_strength": 0.8},
        {"id": "soft_strength_0.60", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 0.6},
        {"id": "soft_strength_1.00", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 1.0},
    )
    specs = build_cpp_feature_map_specs(
        {"z15_26": grid, "z45_56": grid, "z75_86": grid},
        ((3, 3, 3), (5, 5, 3), (7, 7, 3)),
        primary,
        sensitivity,
    )
    protocol = validate_cpp_feature_map_protocol(specs, expected_primary_count=162, expected_sensitivity_count=324)
    assert protocol == {"map_count": 486, "primary_count": 162, "sensitivity_count": 324, "dataset_count": 3, "window_count": 3}


def test_formal_weight_request_rejects_legacy_joint_geometric_mean():
    with pytest.raises(ValueError, match="joint_geometric_mean"):
        _experiment_weight_switches({
            "mode": "soft_high_porosity",
            "enable_soft_high_porosity": True,
            "enable_joint_geometric_mean": True,
        })


def test_formal_weight_request_rejects_unnormalized_cumulant_weights():
    with pytest.raises(ValueError, match="normalized weights"):
        _experiment_weight_switches({
            "mode": "soft_high_porosity",
            "normalize": False,
        })
