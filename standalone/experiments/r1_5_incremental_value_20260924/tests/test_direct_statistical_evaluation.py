from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from cli_runner import CppFeatureMapSpec, build_direct_cpp_feature_map_specs, validate_direct_cpp_feature_map_protocol
from direct_statistical_evaluation import (
    DEFAULT_BASELINES,
    DEFAULT_HIGH_FEATURES,
    DirectCase,
    apply_common_support,
    cliffs_delta,
    compute_direct_effects,
    compute_single_map_advantages,
    derive_standardized_cumulants,
    joint_spatial_block_bootstrap,
    lambda4_equivalence_error,
    standardized_group_contrast,
)
from grid_io import read_grid, write_map
from report_direct_results import _map_slice
from window_support import extract_complete_windows


def _load_direct_orchestrator():
    path = Path(__file__).resolve().parents[2] / "run_direct_sections_44_46.py"
    specification = importlib.util.spec_from_file_location("direct_sections_orchestrator_test", path)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _frame(offset: float = 0.0) -> pd.DataFrame:
    coordinates = [(x, y, 0) for y in range(4) for x in range(4)]
    values = np.arange(len(coordinates), dtype=float) + offset
    channel = np.asarray([(x + y) >= 3 for x, y, _ in coordinates], dtype=bool)
    variance = 1.0 + values / 100.0
    frame = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    frame["channel"] = channel
    frame["channel_occupancy"] = channel.astype(float) * 0.75 + (~channel).astype(float) * 0.25
    frame["analysis_valid"] = True
    frame["mean"] = values
    frame["variance"] = variance
    frame["high_fraction"] = values / values.max()
    frame["lambda3"] = values / 10.0
    frame["lambda4"] = values / 8.0
    frame["kappa4"] = values / 4.0
    frame["lambda5"] = values / 6.0
    frame["lambda6"] = values / 5.0
    return frame


def _case(weighting: str, *, offset: float = 0.0) -> DirectCase:
    return DirectCase("4.6", "z45_56", "5x5x3", weighting, weighting, "primary", _frame(offset))


def _grid(tmp_path: Path):
    coordinates = pd.DataFrame([(x, y, z) for z in range(3) for y in range(7) for x in range(7)], columns=["x", "y", "z"])
    source = tmp_path / "direct_protocol.grid"
    write_map(source, coordinates, 0.1 + 0.01 * coordinates["x"] + 0.02 * coordinates["y"] + 0.03 * coordinates["z"], "value")
    return read_grid(source, "value")


def _direct_specs(tmp_path: Path):
    grid = _grid(tmp_path)
    return build_direct_cpp_feature_map_specs(
        {"z15_26": grid, "z45_56": grid, "z75_86": grid},
        section45_windows=((3, 3, 3), (5, 5, 3), (7, 7, 3)),
        section46_primary_window=(5, 5, 3),
        uniform_weight={"id": "uniform", "mode": "uniform"},
        nonuniform_primary_weights=(
            {"id": "coordinate_gaussian_sigma_1.00", "mode": "coordinate_gaussian", "bandwidth": 1.0},
            {"id": "soft_high_porosity_width_0.02_strength_0.80", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 0.8},
        ),
        sensitivity_weights=(
            {"id": "coordinate_gaussian_sigma_0.75", "mode": "coordinate_gaussian", "bandwidth": 0.75},
            {"id": "coordinate_gaussian_sigma_1.25", "mode": "coordinate_gaussian", "bandwidth": 1.25},
            {"id": "soft_width_0.01", "mode": "soft_high_porosity", "soft_transition_width": 0.01, "soft_weight_strength": 0.8},
            {"id": "soft_width_0.04", "mode": "soft_high_porosity", "soft_transition_width": 0.04, "soft_weight_strength": 0.8},
            {"id": "soft_strength_0.60", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 0.6},
            {"id": "soft_strength_1.00", "mode": "soft_high_porosity", "soft_transition_width": 0.02, "soft_weight_strength": 1.0},
        ),
    )


def test_direct_map_matrix_has_the_frozen_102_unique_identities(tmp_path: Path) -> None:
    specs = _direct_specs(tmp_path)
    protocol = validate_direct_cpp_feature_map_protocol(specs)
    identities = {
        (spec.section, spec.dataset_id, spec.window_size, spec.weighting_id, spec.parameter_set, spec.parameter_role, spec.order)
        for spec in specs
    }
    assert protocol["map_count"] == 102
    assert protocol["primary_count"] == 90
    assert protocol["sensitivity_count"] == 12
    assert len(identities) == 102
    with pytest.raises(ValueError, match="unexpected counts"):
        validate_direct_cpp_feature_map_protocol(specs[:-1])


def test_analysis_smoke_matrix_exercises_primary_and_sidewing_paths(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    specs = orchestrator._analysis_smoke_specs(_direct_specs(tmp_path))
    assert len(specs) == 34
    assert sum(spec.parameter_role == "primary" for spec in specs) == 30
    assert sum(spec.parameter_role == "sensitivity" for spec in specs) == 4
    assert {spec.dataset_id for spec in specs} == {"z45_56"}


def test_frozen_direct_configs_keep_their_logical_source_names(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    sources = {}
    for label in ("template_component_config", "controlled_component_config"):
        source = tmp_path / label / "experiment_config.json"
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(json.dumps({"label": label}), encoding="utf-8")
        sources[label] = source
    copied = orchestrator._copy_frozen_configs(tmp_path / "run", sources)
    assert {path.name for path in copied.values()} == {
        "template_component_config.json", "controlled_component_config.json"
    }
    assert json.loads(copied["template_component_config"].read_text(encoding="utf-8"))["label"] == "template_component_config"
    assert json.loads(copied["controlled_component_config"].read_text(encoding="utf-8"))["label"] == "controlled_component_config"


def test_formal_bootstrap_case_counts_scale_from_the_full_102_map_matrix(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    assert orchestrator._formal_bootstrap_case_counts(_direct_specs(tmp_path)) == {
        "primary": 12,
        "scale": 6,
        "large_block": 12,
        "parameter_sensitivity": 8,
    }


def test_time_projection_includes_case_scaled_bootstrap_and_full_workflow_components(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    manifest = pd.DataFrame(
        {
            "dataset_id": ["z45_56", "z45_56", "z45_56"],
            "window": ["5x5x3", "5x5x3", "5x5x3"],
            "weighting": ["uniform", "uniform", "uniform"],
            "order": [2, 5, 6],
            "cli_elapsed_seconds": [2.0, 2.0, 2.0],
        }
    )
    diagnostics = {
        "primary": {"elapsed_seconds": 4.0, "replicates": 100, "case_count": 4},
        "scale": {"elapsed_seconds": 2.0, "replicates": 100, "case_count": 2},
        "large_block": {"elapsed_seconds": 3.0, "replicates": 100, "case_count": 4},
        "parameter_sensitivity": {"elapsed_seconds": 1.0, "replicates": 100, "case_count": 4},
    }
    common = {
        "bootstrap": {"primary_replicates": 2000, "scale_replicates": 1000, "large_block_replicates": 500},
        "formal_time_budget_seconds": 21600,
    }
    projection = orchestrator._write_time_projection(
        tmp_path,
        manifest,
        mode="analysis-smoke",
        matrix_elapsed_seconds=6.0,
        bootstrap_diagnostics=diagnostics,
        common=common,
        component_seconds={
            "release_build_and_tests": 10.0,
            "preprocessing": 5.0,
            "section44": 2.0,
            "section44_formal_scale": 2.0,
            "controlled_positive_controls": 1.0,
            "controlled_positive_controls_formal_scale": 3.0,
            "cpp_python_parity": 4.0,
            "direct_analysis_overhead": 3.0,
            "report_generation": 2.0,
            "evidence_sealing_estimate": 1.0,
        },
        formal_bootstrap_case_counts={"primary": 12, "scale": 6, "large_block": 12, "parameter_sensitivity": 8},
    )
    assert projection["bootstrap_projection_by_scope"]["primary"]["case_scale_factor"] == pytest.approx(3.0)
    assert projection["bootstrap_projection_by_scope"]["parameter_sensitivity"]["case_scale_factor"] == pytest.approx(2.0)
    assert projection["formal_fixed_component_projection_seconds"]["section44"] == pytest.approx(4.0)
    assert projection["formal_fixed_component_projection_seconds"]["controlled_positive_controls"] == pytest.approx(3.0)
    assert projection["formal_projection"]["estimated_total_seconds"] > (
        projection["formal_projection"]["estimated_map_seconds"]
        + projection["formal_projection"]["estimated_bootstrap_seconds"]
    )
    assert projection["formal_projection"]["time_budget_pass"]


def test_common_support_manifest_keeps_only_explicit_common_rows_and_masks(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    support = pd.DataFrame(
        [
            {"row_type": "raw_support", "support_scope": "individual_weighting_before_common_support", "section": "4.6", "parameter_role": "primary", "weighting": "coordinate_gaussian"},
            {"row_type": "common_support", "support_scope": "section46_all_three_primary_common", "section": "4.6", "parameter_role": "primary", "weighting": "all_three"},
            {"row_type": "support_mask", "support_scope": "section46_all_three_primary_common", "section": "4.6", "parameter_role": "primary", "weighting": "uniform"},
            {"row_type": "support_mask", "support_scope": "section45_individual", "section": "4.5", "parameter_role": "primary", "weighting": "uniform"},
        ]
    )
    support.to_csv(tmp_path / "support_summary.csv", index=False)
    orchestrator._common_support_manifest(tmp_path)
    result = pd.read_csv(tmp_path / "common_support_manifest.csv")
    assert set(result["row_type"]) == {"common_support", "support_mask"}
    assert set(result["support_scope"]) == {"section46_all_three_primary_common"}


def test_map_slice_masks_invalid_boolean_facies_centres() -> None:
    values = pd.DataFrame(
        {
            "x": [0, 1, 0, 1],
            "y": [0, 0, 1, 1],
            "z": [50, 50, 50, 50],
            "channel": [False, True, True, False],
            "analysis_valid": [True, False, True, True],
        }
    )
    image, _, _ = _map_slice(values, "channel", 50)
    assert np.isnan(image[0, 1])
    assert image[1, 0] == pytest.approx(1.0)


def test_independent_raw_lambda4_uses_window_values_and_weights(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    grid = _grid(tmp_path)
    spec = CppFeatureMapSpec(
        dataset_id="z45_56",
        grid=grid,
        window_size=(3, 3, 3),
        order=4,
        weight_config={"id": "uniform", "mode": "uniform"},
        section="4.5",
        weighting_id="uniform",
        parameter_set="uniform",
    )
    coordinates = extract_complete_windows(grid, (3, 3, 3), (1, 1, 1)).coordinates
    raw = orchestrator._independent_raw_lambda4(grid, spec, coordinates, variance_epsilon=1.0e-12)
    values = extract_complete_windows(grid, (3, 3, 3), (1, 1, 1)).values[0]
    centered = values - values.mean()
    expected = float(np.mean(centered**4) / np.mean(centered**2) ** 2 - 3.0)
    assert float(raw.loc[0, "python_raw_excess_kurtosis"]) == pytest.approx(expected)


def test_lambda4_is_explicitly_identical_to_excess_kurtosis() -> None:
    frame = pd.DataFrame(
        {
            "x": [0, 1], "y": [0, 0], "z": [0, 0],
            "kappa1": [0.1, 0.2], "kappa2": [2.0, 4.0], "kappa3": [0.5, -0.5],
            "kappa4": [6.0, -8.0], "kappa5": [1.0, 2.0], "kappa6": [3.0, 4.0],
        }
    )
    result = derive_standardized_cumulants(frame, 1.0e-12)
    audit = lambda4_equivalence_error(result)
    assert np.allclose(result["lambda4"], result["excess_kurtosis"])
    assert audit["finite"]
    assert audit["max_abs_error"] <= 1.0e-12


def test_standardized_group_contrast_and_cliffs_delta_handle_ties() -> None:
    values = np.asarray([1.0, 2.0, 2.0, 3.0])
    channel = np.asarray([False, True, False, True])
    expected = (np.mean([2.0, 3.0]) - np.mean([1.0, 2.0])) / np.sqrt((np.var([2.0, 3.0], ddof=1) + np.var([1.0, 2.0], ddof=1)) / 2.0)
    assert standardized_group_contrast(values, channel) == pytest.approx(expected)
    delta = cliffs_delta(values, channel)
    auc = 0.875  # Three wins and one tie among four channel/non-channel pairs.
    assert delta == pytest.approx(0.75)
    assert delta == pytest.approx(2.0 * auc - 1.0)


def test_joint_bootstrap_uses_cached_statistics_and_recomputes_the_strongest_baseline() -> None:
    case = _case("uniform")
    effects, _ = compute_direct_effects((case,), (*DEFAULT_BASELINES, *DEFAULT_HIGH_FEATURES))
    advantages = compute_single_map_advantages(effects)
    expected = float(effects.loc[effects["feature"] == "lambda5", "separation_S"].iloc[0] - effects.loc[effects["feature"].isin(DEFAULT_BASELINES), "separation_S"].max())
    assert float(advantages.loc[advantages["feature"] == "lambda5", "advantage_A"].iloc[0]) == pytest.approx(expected)
    bootstrap = joint_spatial_block_bootstrap(
        (case,), features=(*DEFAULT_BASELINES, *DEFAULT_HIGH_FEATURES), block_shape_xyz=(20, 20, 20), replicates=20, seed=31
    )
    row = bootstrap.summary.loc[
        (bootstrap.summary["metric_family"] == "direct_effect")
        & (bootstrap.summary["metric_name"] == "effect_E")
        & (bootstrap.summary["feature"] == "mean")
        & (bootstrap.summary["dataset_id"] == "z45_56")
    ].iloc[0]
    direct = float(effects.loc[effects["feature"] == "mean", "effect_E"].iloc[0])
    assert float(row.point_estimate) == pytest.approx(direct)
    assert float(row.bootstrap_mean) == pytest.approx(direct)
    assert bootstrap.diagnostics["accepted_replicates"] == 20


def test_common_support_is_shared_before_weighting_effects_are_bootstrapped() -> None:
    uniform, coordinate, soft = apply_common_support((_case("uniform"), _case("coordinate_gaussian", offset=0.5), _case("soft_high_porosity", offset=1.0)))[0]
    coordinates = uniform.frame.loc[uniform.frame["analysis_valid"], ["x", "y", "z"]]
    for case in (coordinate, soft):
        assert case.frame.loc[case.frame["analysis_valid"], ["x", "y", "z"]].equals(coordinates)
    result = joint_spatial_block_bootstrap(
        (uniform, coordinate, soft), features=(*DEFAULT_BASELINES, *DEFAULT_HIGH_FEATURES), block_shape_xyz=(20, 20, 20), replicates=10, seed=41
    )
    assert set(result.summary.loc[result.summary["metric_family"] == "weighting_effect", "metric_name"]) == {"W", "G"}


def test_direct_output_contract_only_rejects_explicit_legacy_artifacts(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    (tmp_path / "cpp_feature_map_manifest.csv").write_text(
        "section,dataset_id,window,weighting_id,parameter_set,parameter_role,order,algorithm,"
        "cumulant_centering,cumulant_centering_version\n"
        "4.5,fixture,3x3x3,uniform,uniform,primary,1,AVL,none,order_ge_2.v1\n",
        encoding="utf-8",
    )
    for name in ("RUN_REPORT.md", "RESULTS_INDEX.md", "evidence_manifest.json", "analysis_threshold_notes.txt"):
        (tmp_path / name).write_text("fixture\n", encoding="utf-8")
    result = orchestrator._direct_output_contract(tmp_path, mode="smoke", sealed=False)
    assert result["forbidden_artifact_count"] == 0
    section44_table = tmp_path / "section_4_4" / "tables" / "multiplicity_decisions.csv"
    section44_table.parent.mkdir(parents=True)
    section44_table.write_text("section44 surrogate omnibus table\n", encoding="utf-8")
    assert orchestrator._direct_output_contract(tmp_path, mode="smoke", sealed=False)["forbidden_artifact_count"] == 0
    (tmp_path / "oof_predictions.csv").write_text("fixture\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="prohibited model artifacts"):
        orchestrator._direct_output_contract(tmp_path, mode="smoke", sealed=False)


def test_direct_command_runner_keeps_unicode_diagnostics_as_a_log(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    record = orchestrator._run_command(
        [sys.executable, "-c", "print('构建诊断')"],
        tmp_path / "command.log.txt",
        cwd=tmp_path,
    )
    assert record["exit_code"] == 0
    assert "stdout" in (tmp_path / "command.log.txt").read_text(encoding="utf-8")


def test_preprocessing_runtime_config_resolves_relative_fixture_before_copying(tmp_path: Path) -> None:
    orchestrator = _load_direct_orchestrator()
    fixture = tmp_path / "co_location_fixture.json"
    fixture.write_text("{}\n", encoding="utf-8")
    source = tmp_path / "preprocessing_config.json"
    source.write_text('{"co_location_fixture_path": "co_location_fixture.json"}\n', encoding="utf-8")
    runtime = orchestrator._preprocessing_runtime_config(source, tmp_path / "run" / "preprocessing")
    assert runtime["co_location_fixture_path"] == str(fixture.resolve())
    assert runtime["output_root"] == str((tmp_path / "run" / "preprocessing").resolve())
