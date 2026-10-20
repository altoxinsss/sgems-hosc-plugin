from hashlib import sha256
from pathlib import Path
import json
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from formal_runner import (
    _frozen_section44_protocol,
    _frozen_preprocessing_axis_mapping,
    _response_matrix,
    _response_records,
    _smoke_cases,
    _top_level_entrypoint_id,
    _validate_preprocessing_manifest,
    _validate_frozen_geometry,
    _write_component_index,
    run,
)
from directional_validation import sha256_file, write_columnar_grid
from run_template_directional_experiment import CliResult, _write_cli_ini


def _write_fixture_manifest(tmp_path: Path, *, porosity_grid_value: float = 0.25) -> tuple[Path, dict[str, object], list[dict[str, object]]]:
    root = tmp_path / "preprocessing"
    root.mkdir(parents=True)
    expected = {
        "porosity": 0.25,
        "channel_count": 5,
        "channel_occupancy": 0.625,
        "channel_strict_majority": 1,
        "channel_label_valid": 1,
    }
    output_records = []
    paths: dict[str, Path] = {}
    for property_name, value in {**expected, "porosity": porosity_grid_value}.items():
        path = root / f"{property_name}.grid"
        write_columnar_grid(
            path,
            np.asarray([[[float(value)]]]),
            property_name,
            x=[1],
            y=[1],
            z=[1],
        )
        paths[property_name] = path
        output_records.append({
            "path": str(path),
            "property_name": property_name,
            "row_count": 1,
            "sha256": sha256_file(path),
        })
    fixture_path = root / "fixture.json"
    fixture_path.write_text(
        json.dumps({
            "schema_version": "stanford-vie-co-location-fixture.v1",
            "coordinate_order": ["x", "y", "z"],
            "fixtures": [{"coordinate": [1, 1, 1], "expected": expected}],
        }),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": "stanford-vie-preprocessing.v2",
        "transform": "reshape(z,y,x)->mean_or_count_2x2x2->reverse_z_and_y",
        "axis_mapping": _frozen_preprocessing_axis_mapping(),
        "coarse_spacing_xyz_m": [50.0, 50.0, 2.0],
        "co_location_fixture": {
            "path": str(fixture_path),
            "sha256": sha256_file(fixture_path),
            "results": [{"coordinate_xyz": [1, 1, 1], "status": "PASS", "observed": expected}],
        },
        "subvolumes": [{
            "subvolume_id": "z1_1",
            "shape_zyx": [1, 1, 1],
            "coordinate_range": {"x": [1, 1], "y": [1, 1], "z": [1, 1]},
            "outputs": output_records,
        }],
    }
    manifest_path = root / "preprocessing_manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    config = {"preprocessing_manifest": str(manifest_path)}
    datasets = [{
        "id": "stanford_vie_z1_1",
        "attributes": [
            {
                "property_id": "porosity",
                "grid_path": str(paths["porosity"]),
                "property_name": "porosity",
                "field_type": "continuous",
                "surrogate_generator": "continuous_iaaft_swap",
                "connectivity_rule": "high_porosity_excursion",
            },
            {
                "property_id": "channel_occupancy",
                "grid_path": str(paths["channel_occupancy"]),
                "property_name": "channel_occupancy",
                "field_type": "discrete",
                "valid_mask_path": str(paths["channel_label_valid"]),
                "valid_mask_property": "channel_label_valid",
                "surrogate_generator": "discrete_h2_constrained_swap_walk",
                "connectivity_rule": "channel_strict_majority",
            },
        ],
    }]
    return manifest_path, config, datasets


def test_formal_runner_rejects_legacy_sampling_before_any_run(tmp_path):
    config = tmp_path / "legacy.json"
    config.write_text('{"template_sampling": "legacy_search"}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="refuses non-exhaustive"):
        run(config, tmp_path / "output", "smoke")


def test_cli_ini_rejects_legacy_sampling_mode(tmp_path):
    with pytest.raises(ValueError, match="refuses legacy"):
        _write_cli_ini(
            tmp_path / "case.ini", tmp_path / "source.grid", "value", "target", tmp_path / "out",
            nod=[2, 1], nlag=[2, 2], lag=[1, 1], azimuth=[0, 90], dip=[0, 0],
            displayed_dirs=[1, 2], coordinate_dimensions=2, source_spacing=[50.0, 50.0, 2.0],
            direction_vectors=[[1, 0], [0, 1]], template_sampling="legacy_search",
        )


def test_confirmatory_response_matrix_excludes_zero_lag_rows():
    rows = []
    for source_kind, surrogate_id, value in (("original", "", 1.0), ("surrogate", 1, 0.5), ("surrogate", 2, 0.75)):
        for lag_cell, positions in ((0, "0;0"), (1, "1;1")):
            rows.append({
                "dimension": "2d", "algorithm": "AVL Based algorithm", "property_id": "porosity",
                "template_id": "axis_0_90", "dataset_id": "d1", "lag_cell": lag_cell,
                "lag_positions": positions, "source_kind": source_kind, "surrogate_id": surrogate_id,
                "standardized_k4": value + lag_cell,
            })
    matrix = _response_matrix(rows, "porosity", {"axis_0_90"}, 2)
    assert matrix.shape == (3, 1)
    assert np.allclose(matrix[:, 0], [2.0, 1.5, 1.75])


def test_frozen_configuration_covers_the_full_signed_direction_protocol():
    config = json.loads((ROOT / "experiment_config.json").read_text(encoding="utf-8"))
    _validate_frozen_geometry(config)
    smoke_2d, smoke_3d = _smoke_cases(config, config["directional_cases"], config["three_dimensional_cases"])
    assert {case["id"] for case in smoke_2d} == {
        "axis_0_90", "axis_90_180", "diag_45_135", "diag_135_225"
    }
    assert all(case["nlag"] == [2, 2] for case in smoke_2d)
    assert smoke_3d[0]["nlag"] == [2, 2, 2]

    missing_rotation_absolute_tolerance = json.loads(json.dumps(config))
    del missing_rotation_absolute_tolerance["rotation_validation"]["max_absolute_error"]
    with pytest.raises(ValueError, match="absolute and relative rotation tolerances"):
        _validate_frozen_geometry(missing_rotation_absolute_tolerance)


def test_frozen_configuration_requires_exact_three_dimensional_xyz_production_geometry():
    config = json.loads((ROOT / "experiment_config.json").read_text(encoding="utf-8"))
    case = config["three_dimensional_cases"][0]
    assert case["id"] == "orthogonal_xyz"
    assert case["integer_direction_vectors"] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    _validate_frozen_geometry(config)

    wrong_axis = json.loads(json.dumps(config))
    wrong_axis["three_dimensional_cases"][0]["integer_direction_vectors"][2] = [1, 0, 0]
    with pytest.raises(ValueError, match="exact X/Y/Z"):
        _validate_frozen_geometry(wrong_axis)

    non_integer = json.loads(json.dumps(config))
    non_integer["three_dimensional_cases"][0]["integer_direction_vectors"][0][0] = 1.0
    with pytest.raises(ValueError, match="must use JSON integers"):
        _validate_frozen_geometry(non_integer)


def test_preprocessing_manifest_validates_axis_mapping_fixture_hash_and_grid_values(tmp_path):
    manifest_path, config, datasets = _write_fixture_manifest(tmp_path / "good")
    _, _, verified, contract = _validate_preprocessing_manifest(
        config,
        manifest_path.parent,
        datasets,
        [50.0, 50.0, 2.0],
    )
    assert len(verified) == 3
    assert contract["axis_mapping"] == "PASS"
    assert contract["fixture_coordinate_count"] == 1
    assert contract["fixture_grid_file_count"] == 5

    altered_axis = json.loads(manifest_path.read_text(encoding="utf-8"))
    altered_axis["axis_mapping"]["source_reduced_to_output_index"]["y"] = "output_y_index = source_reduced_y_index"
    altered_axis_path = manifest_path.parent / "altered_axis_manifest.json"
    altered_axis_path.write_text(json.dumps(altered_axis), encoding="utf-8")
    altered_axis_config = {"preprocessing_manifest": str(altered_axis_path)}
    with pytest.raises(RuntimeError, match="axis_mapping"):
        _validate_preprocessing_manifest(altered_axis_config, altered_axis_path.parent, datasets, [50.0, 50.0, 2.0])

    bad_manifest_path, bad_config, bad_datasets = _write_fixture_manifest(tmp_path / "bad_grid", porosity_grid_value=0.35)
    with pytest.raises(RuntimeError, match="porosity does not match"):
        _validate_preprocessing_manifest(bad_config, bad_manifest_path.parent, bad_datasets, [50.0, 50.0, 2.0])


def test_formal_component_requires_top_level_entrypoint_binding():
    with pytest.raises(ValueError, match="top_level_entrypoint_id"):
        _top_level_entrypoint_id({}, "formal")
    assert _top_level_entrypoint_id({}, "smoke") == ""
    assert _top_level_entrypoint_id({"top_level_entrypoint_id": "formal-44-46-test"}, "formal") == "formal-44-46-test"


def test_formal_component_rejects_direct_execution_before_any_component_work(tmp_path):
    config = tmp_path / "component_without_top_level_binding.json"
    config.write_text('{"template_sampling": "exhaustive_integer_offsets"}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="top_level_entrypoint_id"):
        run(config, tmp_path / "output", "formal")


def test_formal_component_uses_a_component_scoped_result_index(tmp_path):
    index = _write_component_index(tmp_path)
    assert index.name == "COMPONENT_RESULTS_INDEX.md"
    assert not (tmp_path / "RESULTS_INDEX.md").exists()
    assert "top-level formal entrypoint" in index.read_text(encoding="utf-8")


def test_frozen_section44_protocol_hash_binds_the_executed_geometry(tmp_path):
    runtime_config = json.loads((ROOT / "experiment_config.json").read_text(encoding="utf-8"))
    frozen_source = tmp_path / "formal_section44.json"
    frozen_source.write_bytes((ROOT.parent / "configs" / "formal_section44.json").read_bytes())
    protocol = json.loads(frozen_source.read_text(encoding="utf-8"))
    runtime_config["formal_section44_protocol"] = protocol
    runtime_config["formal_section44_protocol_source_path"] = str(frozen_source)
    runtime_config["formal_section44_protocol_sha256"] = sha256(frozen_source.read_bytes()).hexdigest()

    resolved, binding = _frozen_section44_protocol(runtime_config, "formal")
    assert resolved == protocol
    assert binding == {"source_path": str(frozen_source.resolve()), "sha256": runtime_config["formal_section44_protocol_sha256"]}
    _validate_frozen_geometry(runtime_config, resolved)

    lag_drift = json.loads(json.dumps(runtime_config))
    lag_drift["directional_cases"][0]["nlag"] = [6, 6]
    with pytest.raises(ValueError, match="two-dimensional case axis_0_90 must use frozen NLAG=\\(5,5\\)"):
        _validate_frozen_geometry(lag_drift, resolved)

    frozen_source.write_text(json.dumps({**protocol, "two_dimensional_nlag": 6}), encoding="utf-8")
    with pytest.raises(RuntimeError, match="SHA-256"):
        _frozen_section44_protocol(runtime_config, "formal")


def test_response_records_use_the_run_directory_for_stable_relative_paths(tmp_path):
    run_dir = tmp_path / "run"
    output = run_dir / "cli_outputs" / "response.txt"
    output.parent.mkdir(parents=True)
    output.write_text("placeholder\n", encoding="utf-8")
    result = CliResult("original", "case", output, output.with_suffix(".meta.txt"), output.with_suffix(".support.csv"), output, output, {})
    rows = _response_records(
        dataset_id="d", property_id="porosity", field_type="continuous", condition="original",
        source_kind="original", surrogate_id=None, template_id="axis_0_90", dimension="2d",
        algorithm="AVL Based algorithm", run_dir=run_dir, result=result,
        rows=[{
            "lag_cell": 0, "lag_positions": "0;1", "integer_offsets": "", "physical_offsets": "",
            "valid_translations": 3, "support_count": 3, "raw_c4": 2.0, "grid_x": 0, "grid_y": 1, "grid_z": 0,
        }], global_kappa2=2.0,
    )
    assert rows[0]["output_grid"] == "cli_outputs/response.txt"
    assert rows[0]["zero_lag_baseline"]
    assert not rows[0]["inference_lag_cell"]
