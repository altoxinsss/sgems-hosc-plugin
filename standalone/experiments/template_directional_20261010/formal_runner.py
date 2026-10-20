from __future__ import annotations

"""Formal Section 4.4 integer-offset template experiment runner.

This module is the only Section 4.4 experiment entry point.  It dispatches
the production ``template_cli`` with ``exhaustive_integer_offsets`` and keeps
the reference, surrogate, and evidence paths separate by construction.
"""

import csv
import copy
import json
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from directional_validation import (
    axis_semivariogram_3d,
    connectivity_metrics,
    directional_semivariogram_2d,
    generate_continuous_h2_surrogates,
    generate_discrete_h2_surrogates,
    h2_offsets_from_templates,
    holm_step_down,
    read_columnar_grid,
    sha256_file,
    signed_integer_direction,
    surrogate_diversity_diagnostics,
    symmetric_leave_one_out_test,
    template_response,
    template_response_nd,
    transform_offsets_2d,
    write_columnar_grid,
)
from run_template_directional_experiment import (
    CliResult,
    _matching_explicit_surface,
    _read_volume,
    _relative_l2,
    _resolve,
    _run_cli,
    _run_native_checks,
    _write_cli_ini,
)


ROOT = Path(__file__).resolve().parent


_FIXTURE_PROPERTY_NAMES = (
    "porosity",
    "channel_count",
    "channel_occupancy",
    "channel_strict_majority",
    "channel_label_valid",
)

_FORMAL_3D_CASE_ID = "orthogonal_xyz"
_FORMAL_3D_VECTORS = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def _frozen_preprocessing_axis_mapping() -> dict[str, Any]:
    """Return the immutable Stanford VI-E coarse-grid coordinate contract."""

    return {
        "input_array_order": ["z", "y", "x"],
        "output_array_order": ["z", "y", "x"],
        "fine_shape_zyx": [200, 200, 150],
        "coarse_shape_zyx": [100, 100, 75],
        "reduction": "contiguous_2x2x2_blocks; porosity=mean; facies=channel_count",
        "source_reduced_to_output_index": {
            "z": "output_z_index = 99 - source_reduced_z_index",
            "y": "output_y_index = 99 - source_reduced_y_index",
            "x": "output_x_index = source_reduced_x_index",
        },
        "output_coordinate_to_index": {
            "x": "x_index = x_coordinate - 1",
            "y": "y_index = y_coordinate - 1",
            "z": "z_index = z_coordinate - 1",
        },
    }


_SECTION44_PROTOCOL_KEYS = frozenset({
    "section",
    "runner",
    "sampling_mode",
    "production_algorithm",
    "audit_algorithm",
    "two_dimensional_nlag",
    "three_dimensional_nlag",
    "formal_surrogates_per_dataset_property",
    "smoke_surrogates_per_dataset_property",
    "surrogate_minimum_pair_difference_fraction",
    "surrogate_maximum_duplicate_fraction",
    "continuous_surrogate",
    "discrete_surrogate",
    "h2_max_relative_error",
    "h2_weighted_relative_rmse",
    "rotation_relative_l2",
    "bell_avl_max_abs",
    "bell_avl_relative_l2",
})

_SECTION44_SURROGATE_GENERATORS = {
    "iaaft_then_variogram_swap": ("continuous", "continuous_iaaft_swap"),
    "histogram_preserving_swap_annealing": ("discrete", "discrete_h2_constrained_swap_walk"),
}


def _frozen_section44_protocol(
    config: Mapping[str, Any],
    mode: str,
) -> tuple[dict[str, Any] | None, dict[str, str] | None]:
    """Load and validate the top-level Section 4.4 protocol binding.

    A formal component run is only valid when the protocol object embedded in
    its runtime configuration is byte-bound to the frozen top-level JSON
    snapshot.  Development smoke runs without that binding remain available
    for local diagnostics, but the top-level smoke and formal paths always
    supply it.
    """

    raw_protocol = config.get("formal_section44_protocol")
    if raw_protocol is None:
        if mode == "formal":
            raise ValueError("formal Section 4.4 execution requires a frozen formal_section44 protocol binding")
        return None, None
    if not isinstance(raw_protocol, Mapping):
        raise ValueError("formal_section44_protocol must be a JSON object")
    protocol = dict(raw_protocol)
    missing = sorted(_SECTION44_PROTOCOL_KEYS.difference(protocol))
    unknown = sorted(set(protocol).difference(_SECTION44_PROTOCOL_KEYS))
    if missing or unknown:
        raise ValueError(f"formal_section44 protocol keys are invalid; missing={missing}, unknown={unknown}")

    source_text = config.get("formal_section44_protocol_source_path")
    expected_hash = config.get("formal_section44_protocol_sha256")
    if not isinstance(source_text, str) or not source_text or not isinstance(expected_hash, str) or not expected_hash:
        raise ValueError("frozen Section 4.4 protocol requires source path and SHA-256 binding")
    source_path = Path(source_text).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"frozen Section 4.4 protocol source does not exist: {source_path}")
    actual_hash = sha256_file(source_path)
    if actual_hash != expected_hash:
        raise RuntimeError("frozen Section 4.4 protocol SHA-256 does not match its runtime binding")
    if _load(source_path) != protocol:
        raise RuntimeError("embedded Section 4.4 protocol differs from its frozen source")

    def integer(name: str) -> int:
        value = protocol[name]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"formal Section 4.4 protocol field {name} must be an integer")
        return int(value)

    def positive_number(name: str) -> float:
        value = protocol[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"formal Section 4.4 protocol field {name} must be numeric")
        result = float(value)
        if not np.isfinite(result) or result <= 0.0:
            raise ValueError(f"formal Section 4.4 protocol field {name} must be positive and finite")
        return result

    if str(protocol["section"]) != "4.4":
        raise ValueError("formal Section 4.4 protocol must identify section 4.4")
    if str(protocol["runner"]) != "template_directional_20261010/run_template_directional_experiment.py":
        raise ValueError("formal Section 4.4 protocol names an unsupported component runner")
    if str(protocol["sampling_mode"]) != "exhaustive_integer_offsets":
        raise ValueError("formal Section 4.4 protocol requires exhaustive integer-offset sampling")
    if str(protocol["production_algorithm"]) != "AVL Based algorithm":
        raise ValueError("formal Section 4.4 protocol requires AVL as the production algorithm")
    if str(protocol["audit_algorithm"]) != "Bell algorithm":
        raise ValueError("formal Section 4.4 protocol requires Bell as the audit algorithm")
    if integer("two_dimensional_nlag") < 2 or integer("three_dimensional_nlag") < 2:
        raise ValueError("formal Section 4.4 protocol requires at least two lag positions in each dimension")
    if integer("formal_surrogates_per_dataset_property") != 199:
        raise ValueError("formal Section 4.4 protocol requires 199 accepted surrogates per dataset/property")
    if integer("smoke_surrogates_per_dataset_property") != 30:
        raise ValueError("formal Section 4.4 protocol requires 30 accepted surrogates per smoke dataset/property")
    for name in (
        "surrogate_minimum_pair_difference_fraction",
        "h2_max_relative_error",
        "h2_weighted_relative_rmse",
    ):
        if positive_number(name) > 1.0:
            raise ValueError(f"formal Section 4.4 protocol field {name} must not exceed one")
    duplicate_fraction = protocol["surrogate_maximum_duplicate_fraction"]
    if isinstance(duplicate_fraction, bool) or not isinstance(duplicate_fraction, (int, float)):
        raise ValueError("formal Section 4.4 protocol field surrogate_maximum_duplicate_fraction must be numeric")
    if not 0.0 <= float(duplicate_fraction) < 1.0:
        raise ValueError("formal Section 4.4 protocol field surrogate_maximum_duplicate_fraction must be in [0, 1)")
    positive_number("rotation_relative_l2")
    positive_number("bell_avl_max_abs")
    positive_number("bell_avl_relative_l2")
    if str(protocol["continuous_surrogate"]) not in _SECTION44_SURROGATE_GENERATORS:
        raise ValueError("formal Section 4.4 protocol declares an unsupported continuous surrogate generator")
    if str(protocol["discrete_surrogate"]) not in _SECTION44_SURROGATE_GENERATORS:
        raise ValueError("formal Section 4.4 protocol declares an unsupported discrete surrogate generator")
    if _SECTION44_SURROGATE_GENERATORS[str(protocol["continuous_surrogate"])][0] != "continuous":
        raise ValueError("formal Section 4.4 continuous surrogate declaration is not continuous")
    if _SECTION44_SURROGATE_GENERATORS[str(protocol["discrete_surrogate"])][0] != "discrete":
        raise ValueError("formal Section 4.4 discrete surrogate declaration is not discrete")
    return protocol, {"source_path": str(source_path), "sha256": actual_hash}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("run_%Y%m%d_%H%M%S")


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("experiment configuration must be an object")
    return payload


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    values = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in values:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(values)


def _parse_semicolon_ints(value: str) -> list[int]:
    return [int(part) for part in value.split(";") if part != ""]


def _vectors(case: dict[str, Any], dimensions: int) -> list[list[int]]:
    configured = case.get("integer_direction_vectors")
    if configured is not None:
        if not isinstance(configured, list):
            raise ValueError(f"case '{case['id']}' integer_direction_vectors must be a JSON array")
        vectors = []
        for vector in configured:
            if not isinstance(vector, list):
                raise ValueError(f"case '{case['id']}' integer direction vectors must be JSON arrays")
            if any(isinstance(item, bool) or not isinstance(item, int) for item in vector):
                raise ValueError(f"case '{case['id']}' integer direction vectors must use JSON integers")
            vectors.append(list(vector))
    elif dimensions == 2:
        azimuth = case.get("azimuth")
        if not isinstance(azimuth, list) or any(isinstance(angle, bool) or not isinstance(angle, int) for angle in azimuth):
            raise ValueError(f"case '{case['id']}' azimuth values must be JSON integers")
        vectors = [list(signed_integer_direction(angle)) for angle in azimuth]
    else:
        raise ValueError(f"3D case '{case['id']}' must declare integer_direction_vectors")
    if len(vectors) != len(case["nod"]) or any(len(vector) != dimensions for vector in vectors):
        raise ValueError(f"case '{case['id']}' has invalid integer direction vectors")
    if any(not any(vector) for vector in vectors):
        raise ValueError(f"case '{case['id']}' includes a zero direction vector")
    return vectors


def _case_template_positions_2d(case: dict[str, Any]) -> list[np.ndarray]:
    """Return all frozen variable-position lists in array (y,x) order."""

    vectors = _vectors(case, 2)
    nlag = [int(value) for value in case["nlag"]]
    steps = [int(value) for value in case["lag_cells"]]
    nod = [int(value) for value in case["nod"]]
    positions: list[np.ndarray] = []

    def recurse(index: int, selected: list[int]) -> None:
        if index == len(nlag):
            points: list[tuple[int, int]] = [(0, 0)]
            for direction, lag_position in enumerate(selected):
                dx = vectors[direction][0] * steps[direction] * lag_position
                dy = vectors[direction][1] * steps[direction] * lag_position
                points.extend([(dy, dx)] * nod[direction])
            positions.append(np.asarray(points, dtype=int))
            return
        for position in range(nlag[index]):
            recurse(index + 1, [*selected, position])

    recurse(0, [])
    return positions


def _h2_for_cases(cases: list[dict[str, Any]]) -> np.ndarray:
    return h2_offsets_from_templates(
        template for case in cases for template in _case_template_positions_2d(case)
    )


def _read_sidecar_rows(result: CliResult, property_name: str, displayed_dirs: list[int]) -> list[dict[str, Any]]:
    grid = read_columnar_grid(result.output_grid, property_name)
    x_index = {int(value): index for index, value in enumerate(grid.x)}
    y_index = {int(value): index for index, value in enumerate(grid.y)}
    z_index = {int(value): index for index, value in enumerate(grid.z)}
    output: list[dict[str, Any]] = []
    with result.support_csv.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            positions = _parse_semicolon_ints(row["lag_positions"])
            x = positions[displayed_dirs[0] - 1]
            y = positions[displayed_dirs[1] - 1] if len(displayed_dirs) > 1 else 0
            z = positions[displayed_dirs[2] - 1] if len(displayed_dirs) > 2 else 0
            value = float(grid.values[z_index[z], y_index[y], x_index[x]])
            output.append({
                "lag_cell": int(row["lag_cell"]),
                "lag_positions": row["lag_positions"],
                "integer_offsets": row["integer_offsets"],
                "physical_offsets": row["physical_offsets"],
                "valid_translations": int(row["valid_translations"]),
                "support_count": int(row["support_count"]),
                "raw_c4": value,
                "grid_x": x,
                "grid_y": y,
                "grid_z": z,
            })
    if len(output) != int(grid.values.size):
        raise RuntimeError(f"sidecar/output row mismatch for {result.output_grid}")
    return output


def _field_variance(values: np.ndarray, mask: np.ndarray) -> float:
    selected = np.asarray(values, dtype=float)[np.asarray(mask, dtype=bool)]
    if selected.size < 2:
        raise RuntimeError("field has fewer than two valid values")
    result = float(np.var(selected, ddof=0))
    if not np.isfinite(result) or result <= 0.0:
        raise RuntimeError("field has nonpositive global variance")
    return result


def _read_attribute(config_dir: Path, dataset: dict[str, Any], attribute: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = _resolve(config_dir, attribute["grid_path"])
    values, x, y, z = _read_volume(path, str(attribute["property_name"]))
    valid = np.isfinite(values)
    mask_path_text = attribute.get("valid_mask_path")
    if mask_path_text:
        mask_path = _resolve(config_dir, str(mask_path_text))
        mask_property = str(attribute.get("valid_mask_property", "channel_label_valid"))
        mask_values, mx, my, mz = _read_volume(mask_path, mask_property)
        if not (np.array_equal(x, mx) and np.array_equal(y, my) and np.array_equal(z, mz) and mask_values.shape == values.shape):
            raise RuntimeError(f"valid mask does not align with {attribute['property_id']}")
        valid &= mask_values > 0.5
    if not np.any(valid):
        raise RuntimeError(f"attribute {attribute['property_id']} has no valid values")
    return values, valid, x, y, z


def _attribute_specs(dataset: dict[str, Any]) -> list[dict[str, Any]]:
    attributes = dataset.get("attributes")
    if not isinstance(attributes, list) or not attributes:
        raise ValueError("formal 4.4 datasets must declare a nonempty attributes list")
    required = {
        "property_id", "grid_path", "property_name", "field_type",
        "surrogate_generator", "connectivity_rule",
    }
    for attribute in attributes:
        missing = required.difference(attribute)
        if missing:
            raise ValueError(f"dataset {dataset['id']} attribute is missing {sorted(missing)}")
    return attributes


def _masked_for_cli(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    output = np.asarray(values, dtype=float).copy()
    output[~np.asarray(mask, dtype=bool)] = np.nan
    return output


def _mask_for_connectivity(values: np.ndarray, valid: np.ndarray, rule: str, threshold: float) -> tuple[np.ndarray, str]:
    if rule == "high_porosity_excursion":
        return (np.asarray(values, dtype=float) >= threshold) & valid, "high_porosity_excursion"
    if rule == "channel_strict_majority":
        return (np.asarray(values, dtype=float) > 0.5) & valid, "channel_strict_majority"
    raise ValueError(f"unsupported connectivity rule: {rule}")


def _run_case(
    *,
    executable: Path,
    run_dir: Path,
    source_grid: Path,
    property_name: str,
    target_property: str,
    dataset_id: str,
    property_id: str,
    condition: str,
    case: dict[str, Any],
    dimensions: int,
    source_spacing: list[float],
    algorithm: str,
    template_sampling: str,
    timeout_seconds: int,
) -> tuple[CliResult, list[dict[str, Any]], Path]:
    case_id = str(case["id"])
    route = "bell" if "bell" in algorithm.lower() else "avl"
    output_dir = run_dir / "cli_outputs" / dataset_id / property_id / condition / case_id / route
    ini_path = run_dir / "configs" / f"{dataset_id}__{property_id}__{condition}__{case_id}__{route}.ini"
    # Each run already has a unique directory. Keeping the file stem short is
    # necessary on Windows because the sidecar suffix otherwise crosses the
    # legacy MAX_PATH boundary while the .txt response still appears valid.
    cli_target_property = "response"
    _write_cli_ini(
        ini_path, source_grid, property_name, cli_target_property, output_dir,
        nod=[int(value) for value in case["nod"]], nlag=[int(value) for value in case["nlag"]],
        lag=[int(value) for value in case["lag_cells"]], azimuth=[int(value) for value in case["azimuth"]],
        dip=[int(value) for value in case["dip"]], displayed_dirs=[int(value) for value in case["displayed_dirs"]],
        coordinate_dimensions=dimensions, source_spacing=source_spacing,
        direction_vectors=_vectors(case, dimensions), algorithm=algorithm,
        template_sampling=template_sampling,
    )
    result = _run_cli(
        executable, ini_path, output_dir, cli_target_property, condition, case_id,
        run_dir / "logs" / dataset_id / property_id, timeout_seconds, run_label=route,
    )
    return result, _read_sidecar_rows(result, cli_target_property, [int(item) for item in case["displayed_dirs"]]), ini_path


def _response_records(
    *,
    dataset_id: str,
    property_id: str,
    field_type: str,
    condition: str,
    source_kind: str,
    surrogate_id: int | None,
    template_id: str,
    dimension: str,
    algorithm: str,
    run_dir: Path,
    result: CliResult,
    rows: list[dict[str, Any]],
    global_kappa2: float,
) -> list[dict[str, Any]]:
    result_rows: list[dict[str, Any]] = []
    for row in rows:
        lag_positions = _parse_semicolon_ints(str(row["lag_positions"]))
        result_rows.append({
            "dataset_id": dataset_id,
            "property_id": property_id,
            "field_type": field_type,
            "condition": condition,
            "source_kind": source_kind,
            "surrogate_id": "" if surrogate_id is None else surrogate_id,
            "template_id": template_id,
            "direction_pair": template_id,
            "dimension": dimension,
            "algorithm": algorithm,
            "calculation_path": "production_template_cli",
            "global_kappa2": global_kappa2,
            "standardized_k4": row["raw_c4"] / (global_kappa2 * global_kappa2),
            "zero_lag_baseline": bool(any(position == 0 for position in lag_positions)),
            "inference_lag_cell": bool(dimension == "2d" and all(position > 0 for position in lag_positions)),
            "output_grid": result.output_grid.relative_to(run_dir).as_posix(),
            **row,
        })
    return result_rows


def _flatten_audit(rows: list[dict[str, object]], dataset_id: str, property_id: str) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        h2_row_keys = ("h2_rows", "occupancy_h2_rows", "label_h2_rows")
        copy = {key: value for key, value in row.items() if key not in h2_row_keys}
        copy.update({
            "dataset_id": dataset_id,
            "property_id": property_id,
            "h2_rows_json": json.dumps(row.get("h2_rows", []), sort_keys=True),
            "occupancy_h2_rows_json": json.dumps(row.get("occupancy_h2_rows", []), sort_keys=True),
            "label_h2_rows_json": json.dumps(row.get("label_h2_rows", []), sort_keys=True),
        })
        output.append(copy)
    return output


def _flatten_h2_diagnostics(rows: list[dict[str, object]], dataset_id: str, property_id: str) -> list[dict[str, Any]]:
    """Expose every accepted or rejected H2 offset outside a JSON cell."""

    output: list[dict[str, Any]] = []
    paths = (
        ("continuous_or_generic", "h2_rows"),
        ("occupancy", "occupancy_h2_rows"),
        ("strict_majority", "label_h2_rows"),
    )
    for audit in rows:
        for statistic_path, key in paths:
            for row in audit.get(key, []):
                flattened = dict(row)
                flattened.update({
                    "dataset_id": dataset_id,
                    "property_id": property_id,
                    "surrogate_id": audit.get("surrogate_id", ""),
                    "attempt": audit.get("attempt", ""),
                    "accepted": audit.get("accepted", False),
                    "generator_type": audit.get("generator_type", ""),
                    "statistic_path": statistic_path,
                })
                output.append(flattened)
    return output


def _response_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["dimension"] != "2d" or row["algorithm"] != "AVL Based algorithm":
            continue
        key = (row["dataset_id"], row["property_id"], row["template_id"], row["direction_pair"], row["lag_cell"])
        grouped[key].append(row)
    output: list[dict[str, Any]] = []
    for key, values in sorted(grouped.items()):
        original = [float(row["standardized_k4"]) for row in values if row["source_kind"] == "original"]
        surrogate = [float(row["standardized_k4"]) for row in values if row["source_kind"] == "surrogate"]
        if len(original) != 1:
            continue
        output.append({
            "dataset_id": key[0], "property_id": key[1], "template_id": key[2], "direction_pair": key[3], "lag_cell": key[4],
            "original_standardized_k4": original[0], "surrogate_count": len(surrogate),
            "surrogate_mean": float(np.mean(surrogate)) if surrogate else float("nan"),
            "surrogate_standard_deviation": float(np.std(surrogate, ddof=0)) if surrogate else float("nan"),
            "surrogate_q025": float(np.quantile(surrogate, 0.025)) if surrogate else float("nan"),
            "surrogate_q975": float(np.quantile(surrogate, 0.975)) if surrogate else float("nan"),
        })
    return output


def _response_matrix(
    rows: list[dict[str, Any]],
    property_id: str,
    cases: set[str],
    count: int,
) -> np.ndarray:
    relevant = [
        row for row in rows
        if row["dimension"] == "2d" and row["algorithm"] == "AVL Based algorithm"
        and row["property_id"] == property_id and row["template_id"] in cases
        and row["source_kind"] in {"original", "surrogate"}
        # Lag zero is retained in the response facts as a geometry audit.  It
        # is not part of the frozen directional surrogate inference.
        and all(position > 0 for position in _parse_semicolon_ints(str(row["lag_positions"])))
    ]
    keys = sorted({
        (row["dataset_id"], row["template_id"], int(row["lag_cell"]))
        for row in relevant
    })
    values: dict[tuple[str, int, tuple[Any, ...]], float] = {}
    for row in relevant:
        candidate = "original" if row["source_kind"] == "original" else f"surrogate_{int(row['surrogate_id']):03d}"
        values[(candidate, 0, (row["dataset_id"], row["template_id"], int(row["lag_cell"])))] = float(row["standardized_k4"])
    matrix: list[list[float]] = []
    for candidate in ["original", *[f"surrogate_{index:03d}" for index in range(1, count + 1)]]:
        vector: list[float] = []
        for key in keys:
            value_key = (candidate, 0, key)
            if value_key not in values:
                raise RuntimeError(f"missing response for {property_id}/{candidate}/{key}")
            vector.append(values[value_key])
        matrix.append(vector)
    return np.asarray(matrix, dtype=float)


def _connectivity_test_matrix(connectivity_rows: list[dict[str, Any]], property_id: str, count: int) -> np.ndarray:
    relevant = [row for row in connectivity_rows if row["property_id"] == property_id and row["dimension"] == "2d"]
    datasets = sorted({str(row["dataset_id"]) for row in relevant})
    values: dict[tuple[str, str], float] = {}
    for row in relevant:
        if row["source_kind"] not in {"original", "surrogate"}:
            continue
        candidate = "original" if row["source_kind"] == "original" else f"surrogate_{int(row['surrogate_id']):03d}"
        values[(candidate, str(row["dataset_id"]))] = float(row["largest_component_fraction_of_active"])
    result = []
    for candidate in ["original", *[f"surrogate_{index:03d}" for index in range(1, count + 1)]]:
        result.append([values[(candidate, dataset)] for dataset in datasets])
    return np.asarray(result, dtype=float)


def _render_figures_from_tables(run_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = run_dir / "figures"
    figures.mkdir(exist_ok=True)
    paths: list[Path] = []
    with (run_dir / "tables" / "template_2d_response_long.csv").open("r", encoding="utf-8", newline="") as handle:
        responses = list(csv.DictReader(handle))
    original = [row for row in responses if row["source_kind"] == "original" and row["algorithm"] == "AVL Based algorithm"]
    if original:
        figure, axis = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
        for (dataset, prop, case), group in _group_rows(original, ("dataset_id", "property_id", "template_id")).items():
            ordered = sorted(group, key=lambda row: int(row["lag_cell"]))
            axis.plot([int(row["lag_cell"]) for row in ordered], [float(row["standardized_k4"]) for row in ordered], marker="o", label=f"{dataset}/{prop}/{case}")
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set_xlabel("Lag cell")
        axis.set_ylabel("Standardized fourth cumulant K4")
        axis.legend(fontsize=6, ncol=2)
        path = figures / "directional_standardized_k4.png"
        figure.savefig(path, dpi=320)
        plt.close(figure)
        paths.append(path)
    with (run_dir / "tables" / "surrogate_response_summary.csv").open("r", encoding="utf-8", newline="") as handle:
        summary = list(csv.DictReader(handle))
    if summary:
        chosen = summary[: min(10, len(summary))]
        figure, axis = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
        x = np.arange(len(chosen))
        mean = np.asarray([float(row["surrogate_mean"]) for row in chosen])
        lower = np.asarray([float(row["surrogate_q025"]) for row in chosen])
        upper = np.asarray([float(row["surrogate_q975"]) for row in chosen])
        original_values = np.asarray([float(row["original_standardized_k4"]) for row in chosen])
        axis.errorbar(x, mean, yerr=np.vstack((mean - lower, upper - mean)), fmt="o", color="#5b6770", label="surrogate 95% interval")
        axis.scatter(x, original_values, color="#b62f28", marker="D", label="original")
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set_xlabel("Grouped directional lag cell")
        axis.set_ylabel("Standardized fourth cumulant K4")
        axis.legend(fontsize=8)
        path = figures / "surrogate_standardized_k4_comparison.png"
        figure.savefig(path, dpi=320)
        plt.close(figure)
        paths.append(path)
    with (run_dir / "tables" / "template_3d_response_long.csv").open("r", encoding="utf-8", newline="") as handle:
        three_d = list(csv.DictReader(handle))
    original_3d = [row for row in three_d if row["source_kind"] == "original" and row["algorithm"] == "AVL Based algorithm"]
    if original_3d:
        figure, axis = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
        values = [float(row["standardized_k4"]) for row in original_3d]
        image = axis.scatter(
            [int(row["grid_x"]) for row in original_3d], [int(row["grid_y"]) for row in original_3d],
            c=values, s=48, cmap="coolwarm", vmin=-max(abs(min(values)), abs(max(values))), vmax=max(abs(min(values)), abs(max(values))),
        )
        axis.set_xlabel("Lag index direction 1")
        axis.set_ylabel("Lag index direction 2")
        figure.colorbar(image, ax=axis, label="K4")
        path = figures / "three_dimensional_standardized_k4.png"
        figure.savefig(path, dpi=320)
        plt.close(figure)
        paths.append(path)
    return paths


def _group_rows(rows: Iterable[dict[str, Any]], fields: tuple[str, ...]) -> dict[tuple[str, ...], list[dict[str, Any]]]:
    result: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        result[tuple(str(row[field]) for field in fields)].append(row)
    return result


def _hash_table(root: Path) -> list[dict[str, str]]:
    return [
        {"path": item.relative_to(root).as_posix(), "sha256": sha256_file(item)}
        for item in sorted(root.rglob("*"))
        if item.is_file() and item.name not in {"evidence_manifest.json", "sha256_manifest.csv"}
    ]


def _write_report(run_dir: Path, status: str, gates: dict[str, bool], surrogate_count: int) -> Path:
    path = run_dir / "RUN_REPORT.md"
    lines = ["# Section 4.4 Directional Template Run", "", f"Status: **{status}**", "", "## Gates", ""]
    lines.extend(f"- `{name}`: `{'PASS' if value else 'FAIL'}`" for name, value in gates.items())
    lines.extend([
        "", "## Frozen Design", "",
        "- Production templates use signed integer offsets and exhaustive valid translations.",
        "- Original fields are audited with Bell and AVL on identical template supports; surrogate fields use AVL only.",
        "- The H2 surrogate gate reads only histogram and semivariogram quantities, not cumulant or connectivity results.",
        f"- Requested accepted surrogates per two-dimensional field: `{surrogate_count}`.",
        "", "## Output", "",
        "- `tables/template_2d_response_long.csv` and `tables/template_3d_response_long.csv` are the production response facts.",
        "- `tables/H2_offsets.csv`, `tables/directional_variogram_long.csv`, `tables/surrogate_acceptance.csv`, and `tables/surrogate_h2_diagnostics.csv` record the null-model gate.",
        "- `tables/bell_avl_accuracy.csv`, `tables/rotation_validation_long.csv`, and `tables/connectivity_metrics.csv` record audit results.",
        "- `evidence_manifest.json` and `sha256_manifest.csv` provide file-level provenance.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_component_index(run_dir: Path) -> Path:
    path = run_dir / "COMPONENT_RESULTS_INDEX.md"
    path.write_text(
        "# Section 4.4 Component Result Placement\n\n"
        "This component index is consumed by the top-level formal entrypoint. "
        "The study-level result index and final PASS status are written only at the top level.\n\n"
        "## 4.4.1 Directional integer-offset responses\n\n"
        "Use `tables/template_2d_response_long.csv`, `tables/template_3d_response_long.csv`, and "
        "`tables/template_definitions.csv`. The response rows retain every lag cell, signed integer offset, "
        "physical offset, and complete-translation support count.\n\n"
        "## 4.4.2 Second-order-matched surrogate comparison\n\n"
        "Use `tables/H2_offsets.csv`, `tables/directional_variogram_long.csv`, "
        "`tables/surrogate_acceptance.csv`, `tables/surrogate_h2_diagnostics.csv`, `tables/surrogate_response_summary.csv`, and "
        "`tables/multiplicity_decisions.csv`. The porosity omnibus is the primary response test; occupancy and "
        "strict-majority connectivity are separate Holm-adjusted secondary tests.\n\n"
        "## 4.4.3 Rotation, algorithm audit, and connectivity\n\n"
        "Use `tables/rotation_validation_long.csv`, `tables/bell_avl_accuracy.csv`, and "
        "`tables/connectivity_metrics.csv`. Explicit line/rectangle/L responses are marked exploratory and do not "
        "enter the omnibus surrogate tests.\n\n"
        "## Reproducibility\n\n"
        "Use generated `configs/`, `inputs/`, `logs/`, `sha256_manifest.csv`, and `evidence_manifest.json`.\n",
        encoding="utf-8",
    )
    return path


def _fixture_float(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise RuntimeError(f"{context} must be a finite number")
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError(f"{context} must be a finite number")
    return result


def _fixture_integer(value: Any, context: str) -> int:
    numeric = _fixture_float(value, context)
    if not numeric.is_integer():
        raise RuntimeError(f"{context} must be an integer")
    return int(numeric)


def _fixture_coordinate(value: Any, context: str) -> tuple[int, int, int]:
    if not isinstance(value, list) or len(value) != 3:
        raise RuntimeError(f"{context} must be [x, y, z]")
    return tuple(_fixture_integer(component, f"{context}[{index}]") for index, component in enumerate(value))


def _validate_fixture_value_mapping(
    expected: Mapping[str, Any],
    observed: Mapping[str, Any],
    *,
    context: str,
) -> None:
    if set(expected) != set(_FIXTURE_PROPERTY_NAMES):
        raise RuntimeError(f"{context} expected values must enumerate every co-located property")
    if set(observed) != set(_FIXTURE_PROPERTY_NAMES):
        raise RuntimeError(f"{context} observed values must enumerate every co-located property")
    if not np.isclose(
        _fixture_float(observed["porosity"], f"{context} observed porosity"),
        _fixture_float(expected["porosity"], f"{context} expected porosity"),
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise RuntimeError(f"{context} porosity does not match the fixed co-location fixture")
    for property_name in ("channel_count", "channel_strict_majority", "channel_label_valid"):
        if _fixture_integer(observed[property_name], f"{context} observed {property_name}") != _fixture_integer(
            expected[property_name], f"{context} expected {property_name}"
        ):
            raise RuntimeError(f"{context} {property_name} does not match the fixed co-location fixture")
    if not np.isclose(
        _fixture_float(observed["channel_occupancy"], f"{context} observed channel_occupancy"),
        _fixture_float(expected["channel_occupancy"], f"{context} expected channel_occupancy"),
        rtol=0.0,
        atol=1.0e-15,
    ):
        raise RuntimeError(f"{context} channel_occupancy does not match the fixed co-location fixture")

    count = _fixture_integer(expected["channel_count"], f"{context} expected channel_count")
    occupancy = _fixture_float(expected["channel_occupancy"], f"{context} expected channel_occupancy")
    strict_majority = _fixture_integer(expected["channel_strict_majority"], f"{context} expected channel_strict_majority")
    label_valid = _fixture_integer(expected["channel_label_valid"], f"{context} expected channel_label_valid")
    if not np.isclose(occupancy, count / 8.0, rtol=0.0, atol=1.0e-15):
        raise RuntimeError(f"{context} channel_occupancy is inconsistent with channel_count")
    if strict_majority != int(count > 4) or label_valid != int(count != 4):
        raise RuntimeError(f"{context} strict-majority fields are inconsistent with channel_count")


def _validate_preprocessing_axis_mapping_and_fixture(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Verify the frozen Y/Z reversal contract and every fixed co-location value.

    The manifest hash protects the declared preprocessing run.  This additional
    check makes the Section 4.4 component reject a manifest whose coordinate
    description or fixture evidence has drifted before template calculations
    begin, and independently reads the registered output grids at each fixture
    coordinate.
    """

    if manifest.get("transform") != "reshape(z,y,x)->mean_or_count_2x2x2->reverse_z_and_y":
        raise RuntimeError("preprocessing transform does not match the frozen Section 4.4 contract")
    if manifest.get("axis_mapping") != _frozen_preprocessing_axis_mapping():
        raise RuntimeError("preprocessing axis_mapping does not match the frozen Section 4.4 coordinate contract")

    fixture_metadata = manifest.get("co_location_fixture")
    if not isinstance(fixture_metadata, Mapping):
        raise RuntimeError("preprocessing manifest lacks co-location fixture metadata")
    fixture_path_text = fixture_metadata.get("path")
    fixture_hash = fixture_metadata.get("sha256")
    if not isinstance(fixture_path_text, str) or not fixture_path_text or not isinstance(fixture_hash, str) or not fixture_hash:
        raise RuntimeError("preprocessing co-location fixture metadata is incomplete")
    fixture_path = Path(fixture_path_text).resolve()
    if not fixture_path.is_file():
        raise FileNotFoundError(f"preprocessing co-location fixture does not exist: {fixture_path}")
    if sha256_file(fixture_path) != fixture_hash:
        raise RuntimeError("preprocessing co-location fixture SHA-256 does not match the manifest")
    fixture_payload = _load(fixture_path)
    if fixture_payload.get("schema_version") != "stanford-vie-co-location-fixture.v1":
        raise RuntimeError("preprocessing co-location fixture schema is not supported")
    if fixture_payload.get("coordinate_order") != ["x", "y", "z"]:
        raise RuntimeError("preprocessing co-location fixture must use x/y/z coordinates")
    fixture_rows = fixture_payload.get("fixtures")
    manifest_rows = fixture_metadata.get("results")
    if not isinstance(fixture_rows, list) or not fixture_rows:
        raise RuntimeError("preprocessing co-location fixture must contain at least one fixed coordinate")
    if not isinstance(manifest_rows, list) or len(manifest_rows) != len(fixture_rows):
        raise RuntimeError("preprocessing manifest co-location results do not cover the fixture coordinates")

    expected_by_coordinate: dict[tuple[int, int, int], dict[str, Any]] = {}
    for index, row in enumerate(fixture_rows):
        if not isinstance(row, Mapping):
            raise RuntimeError("preprocessing co-location fixture rows must be JSON objects")
        coordinate = _fixture_coordinate(row.get("coordinate"), f"co-location fixture row {index}")
        expected = row.get("expected")
        if coordinate in expected_by_coordinate or not isinstance(expected, Mapping):
            raise RuntimeError("preprocessing co-location fixture coordinates must be unique and have expected values")
        _validate_fixture_value_mapping(expected, expected, context=f"co-location fixture {coordinate}")
        expected_by_coordinate[coordinate] = dict(expected)

    observed_by_coordinate: dict[tuple[int, int, int], dict[str, Any]] = {}
    for index, row in enumerate(manifest_rows):
        if not isinstance(row, Mapping):
            raise RuntimeError("preprocessing co-location result rows must be JSON objects")
        coordinate = _fixture_coordinate(row.get("coordinate_xyz"), f"co-location manifest result {index}")
        observed = row.get("observed")
        if str(row.get("status")) != "PASS" or coordinate in observed_by_coordinate or not isinstance(observed, Mapping):
            raise RuntimeError("preprocessing co-location results must be unique PASS records with observed values")
        observed_by_coordinate[coordinate] = dict(observed)
    if set(observed_by_coordinate) != set(expected_by_coordinate):
        raise RuntimeError("preprocessing co-location result coordinates differ from the frozen fixture")
    for coordinate, expected in expected_by_coordinate.items():
        _validate_fixture_value_mapping(
            expected,
            observed_by_coordinate[coordinate],
            context=f"preprocessing co-location manifest result {coordinate}",
        )

    output_records: dict[tuple[str, str], Mapping[str, Any]] = {}
    subvolume_ranges: list[tuple[str, dict[str, tuple[int, int]], tuple[int, int, int]]] = []
    for subvolume in manifest.get("subvolumes", []):
        if not isinstance(subvolume, Mapping):
            raise RuntimeError("preprocessing manifest subvolume metadata must be JSON objects")
        subvolume_id = str(subvolume.get("subvolume_id", ""))
        coordinate_range = subvolume.get("coordinate_range")
        shape = subvolume.get("shape_zyx")
        if not subvolume_id or not isinstance(coordinate_range, Mapping) or not isinstance(shape, list) or len(shape) != 3:
            raise RuntimeError("preprocessing manifest subvolume metadata is incomplete")
        ranges: dict[str, tuple[int, int]] = {}
        for axis in ("x", "y", "z"):
            values = coordinate_range.get(axis)
            if not isinstance(values, list) or len(values) != 2:
                raise RuntimeError(f"preprocessing subvolume {subvolume_id} has an invalid {axis} coordinate range")
            start = _fixture_integer(values[0], f"preprocessing subvolume {subvolume_id} {axis} range start")
            end = _fixture_integer(values[1], f"preprocessing subvolume {subvolume_id} {axis} range end")
            if end < start:
                raise RuntimeError(f"preprocessing subvolume {subvolume_id} has an inverted {axis} coordinate range")
            ranges[axis] = (start, end)
        expected_shape = (ranges["z"][1] - ranges["z"][0] + 1, ranges["y"][1] - ranges["y"][0] + 1, ranges["x"][1] - ranges["x"][0] + 1)
        actual_shape = tuple(_fixture_integer(value, f"preprocessing subvolume {subvolume_id} shape") for value in shape)
        if actual_shape != expected_shape:
            raise RuntimeError(f"preprocessing subvolume {subvolume_id} shape does not match its coordinate ranges")
        subvolume_ranges.append((subvolume_id, ranges, actual_shape))
        outputs = subvolume.get("outputs")
        if not isinstance(outputs, list):
            raise RuntimeError(f"preprocessing subvolume {subvolume_id} lacks output metadata")
        for output in outputs:
            if not isinstance(output, Mapping):
                raise RuntimeError(f"preprocessing subvolume {subvolume_id} has an invalid output record")
            property_name = str(output.get("property_name", ""))
            key = (subvolume_id, property_name)
            if property_name not in _FIXTURE_PROPERTY_NAMES or key in output_records:
                raise RuntimeError(f"preprocessing subvolume {subvolume_id} has invalid or duplicate property output metadata")
            output_records[key] = output

    grid_cache: dict[tuple[str, str], Any] = {}
    for coordinate, expected in expected_by_coordinate.items():
        x, y, z = coordinate
        matching_subvolumes = [
            (subvolume_id, ranges, shape)
            for subvolume_id, ranges, shape in subvolume_ranges
            if ranges["x"][0] <= x <= ranges["x"][1]
            and ranges["y"][0] <= y <= ranges["y"][1]
            and ranges["z"][0] <= z <= ranges["z"][1]
        ]
        if len(matching_subvolumes) != 1:
            raise RuntimeError(f"co-location fixture {coordinate} must belong to exactly one preprocessing subvolume")
        subvolume_id, ranges, expected_shape = matching_subvolumes[0]
        actual: dict[str, Any] = {}
        for property_name in _FIXTURE_PROPERTY_NAMES:
            key = (subvolume_id, property_name)
            output = output_records.get(key)
            if output is None:
                raise RuntimeError(f"preprocessing subvolume {subvolume_id} lacks {property_name} fixture output metadata")
            if key not in grid_cache:
                output_path = Path(str(output.get("path", ""))).resolve()
                expected_hash = str(output.get("sha256", ""))
                if not output_path.is_file() or not expected_hash or sha256_file(output_path) != expected_hash:
                    raise RuntimeError(f"preprocessing fixture output hash does not match manifest: {output_path}")
                grid = read_columnar_grid(output_path, property_name)
                if tuple(int(value) for value in grid.values.shape) != expected_shape:
                    raise RuntimeError(f"preprocessing fixture output shape differs from its manifest: {output_path}")
                expected_axes = {
                    "x": np.arange(ranges["x"][0], ranges["x"][1] + 1, dtype=int),
                    "y": np.arange(ranges["y"][0], ranges["y"][1] + 1, dtype=int),
                    "z": np.arange(ranges["z"][0], ranges["z"][1] + 1, dtype=int),
                }
                if not (
                    np.array_equal(grid.x, expected_axes["x"])
                    and np.array_equal(grid.y, expected_axes["y"])
                    and np.array_equal(grid.z, expected_axes["z"])
                    and int(grid.values.size) == _fixture_integer(output.get("row_count"), f"preprocessing output {output_path} row_count")
                ):
                    raise RuntimeError(f"preprocessing fixture output coordinates differ from its manifest: {output_path}")
                grid_cache[key] = grid
            grid = grid_cache[key]
            actual[property_name] = float(
                grid.values[
                    int(z - ranges["z"][0]),
                    int(y - ranges["y"][0]),
                    int(x - ranges["x"][0]),
                ]
            )
        _validate_fixture_value_mapping(expected, actual, context=f"preprocessing fixture grid values {coordinate}")

    return {
        "axis_mapping": "PASS",
        "fixture_path": str(fixture_path),
        "fixture_sha256": str(fixture_hash),
        "fixture_coordinate_count": len(expected_by_coordinate),
        "fixture_grid_file_count": len(grid_cache),
    }


def _validate_preprocessing_manifest(
    config: dict[str, Any],
    config_dir: Path,
    datasets: list[dict[str, Any]],
    source_spacing: list[float],
) -> tuple[Path, dict[str, Any], list[dict[str, str]], dict[str, Any]]:
    """Confirm that every configured Section 4.4 input is in the frozen data manifest."""

    manifest_path = _resolve(config_dir, str(config["preprocessing_manifest"]))
    if not manifest_path.exists():
        raise FileNotFoundError(f"preprocessing manifest does not exist: {manifest_path}")
    manifest = _load(manifest_path)
    if str(manifest.get("schema_version")) != "stanford-vie-preprocessing.v2":
        raise RuntimeError("Section 4.4 requires the v2 Stanford VI-E preprocessing manifest")
    manifest_spacing = [float(value) for value in manifest.get("coarse_spacing_xyz_m", [])]
    if manifest_spacing != source_spacing:
        raise RuntimeError("preprocessing spacing does not match the frozen Section 4.4 configuration")
    preprocessing_contract = _validate_preprocessing_axis_mapping_and_fixture(manifest)

    output_hashes: dict[str, str] = {}
    subvolume_ids: set[str] = set()
    for subvolume in manifest.get("subvolumes", []):
        subvolume_ids.add(str(subvolume.get("subvolume_id", "")))
        for output in subvolume.get("outputs", []):
            output_hashes[str(Path(str(output["path"])).resolve())] = str(output["sha256"])
    verified: list[dict[str, str]] = []
    for dataset in datasets:
        dataset_id = str(dataset["id"])
        expected_subvolume = dataset_id.removeprefix("stanford_vie_")
        if expected_subvolume not in subvolume_ids:
            raise RuntimeError(f"dataset {dataset_id} is absent from the preprocessing manifest")
        for attribute in _attribute_specs(dataset):
            entries = [("grid", _resolve(config_dir, str(attribute["grid_path"])))]
            if attribute.get("valid_mask_path"):
                entries.append(("valid_mask", _resolve(config_dir, str(attribute["valid_mask_path"]))))
            for role, path in entries:
                resolved = str(path.resolve())
                if resolved not in output_hashes:
                    raise RuntimeError(f"{role} input is not registered in preprocessing manifest: {path}")
                actual = sha256_file(path)
                if actual != output_hashes[resolved]:
                    raise RuntimeError(f"{role} input hash differs from preprocessing manifest: {path}")
                verified.append({"dataset_id": dataset_id, "property_id": str(attribute["property_id"]), "role": role,
                                 "path": resolved, "sha256": actual})
    return manifest_path, manifest, verified, preprocessing_contract


def _smoke_cases(config: dict[str, Any], cases_2d: list[dict[str, Any]], cases_3d: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Create a bounded but geometrically representative smoke matrix.

    The copied cases keep the frozen signed-direction definitions while using
    two lag positions per axis.  The formal branch always uses the untouched
    configuration, including NLAG=5 and NLAG=4.
    """

    required = set(str(value) for value in config["confirmatory_direction_pairs"])
    for pair in config["rotation_validation"]["pairs"]:
        if bool(pair.get("include_in_smoke", False)):
            required.add(str(pair["original_case_id"]))
            required.add(str(pair["rotated_case_id"]))
    selected_2d = [copy.deepcopy(case) for case in cases_2d if str(case["id"]) in required]
    if not selected_2d:
        raise RuntimeError("smoke matrix selected no two-dimensional directional cases")
    for case in selected_2d:
        case["nlag"] = [min(2, int(value)) for value in case["nlag"]]
    selected_3d = [copy.deepcopy(case) for case in cases_3d[:1]]
    for case in selected_3d:
        case["nlag"] = [min(2, int(value)) for value in case["nlag"]]
    return selected_2d, selected_3d


def _accepted_surrogate_keys(rows: list[dict[str, Any]]) -> set[tuple[str, str, int]]:
    return {
        (str(row["dataset_id"]), str(row["property_id"]), int(row["surrogate_id"]))
        for row in rows
        if bool(row.get("accepted"))
    }


def _validate_frozen_geometry(
    config: dict[str, Any],
    protocol: Mapping[str, Any] | None = None,
) -> None:
    """Reject a configuration that silently drifts from the Section 4.4 protocol."""

    sampling_mode = str(protocol["sampling_mode"]) if protocol is not None else "exhaustive_integer_offsets"
    if str(config.get("template_sampling")) != sampling_mode:
        raise ValueError("formal Section 4.4 runner requires exhaustive integer-offset sampling")
    surrogate_spec = dict(config.get("surrogates", {}))
    if int(surrogate_spec.get("swap_iterations", 0)) <= 0:
        raise ValueError("formal Section 4.4 runner requires positive discrete swap iterations")
    minimum_original_difference = float(surrogate_spec.get("minimum_original_difference_fraction", 0.0))
    if not 0.0 < minimum_original_difference <= 1.0:
        raise ValueError("formal Section 4.4 runner requires a frozen original-to-surrogate difference threshold in (0, 1]")
    minimum_pair_difference = float(surrogate_spec.get("minimum_pair_difference_fraction", 0.0))
    maximum_duplicate_fraction = float(surrogate_spec.get("maximum_duplicate_fraction", 1.0))
    if not 0.0 < minimum_pair_difference <= 1.0:
        raise ValueError("formal Section 4.4 runner requires a frozen positive pairwise-surrogate diversity threshold")
    if not 0.0 <= maximum_duplicate_fraction < 1.0:
        raise ValueError("formal Section 4.4 runner requires a frozen duplicate-fraction threshold in [0, 1)")
    rotation = dict(config.get("rotation_validation", {}))
    if float(rotation.get("max_absolute_error", -1.0)) <= 0.0 or float(rotation.get("max_relative_l2_error", -1.0)) <= 0.0:
        raise ValueError("formal Section 4.4 runner requires positive absolute and relative rotation tolerances")
    two_d = list(config.get("directional_cases", []))
    three_d = list(config.get("three_dimensional_cases", []))
    if not two_d or not three_d:
        raise ValueError("the frozen protocol requires two- and three-dimensional template cases")
    expected_directions = {
        (1, 0), (1, 1), (0, 1), (-1, 1),
        (-1, 0), (-1, -1), (0, -1), (1, -1),
    }
    seen_directions: set[tuple[int, int]] = set()
    case_ids: set[str] = set()
    two_dimensional_nlag = int(protocol["two_dimensional_nlag"]) if protocol is not None else 5
    three_dimensional_nlag = int(protocol["three_dimensional_nlag"]) if protocol is not None else 4
    for case in two_d:
        case_id = str(case["id"])
        if case_id in case_ids:
            raise ValueError(f"duplicate directional case id: {case_id}")
        case_ids.add(case_id)
        if [int(value) for value in case["nlag"]] != [two_dimensional_nlag, two_dimensional_nlag] or [int(value) for value in case["lag_cells"]] != [1, 1]:
            raise ValueError(
                f"two-dimensional case {case_id} must use frozen NLAG=({two_dimensional_nlag},{two_dimensional_nlag}) and LAG=(1,1)"
            )
        vectors = _vectors(case, 2)
        azimuth = case.get("azimuth")
        if not isinstance(azimuth, list) or len(azimuth) != len(vectors):
            raise ValueError(f"two-dimensional case {case_id} must declare one azimuth per integer direction")
        expected_vectors = [
            signed_integer_direction(angle)
            for angle in azimuth
        ]
        if tuple(tuple(vector) for vector in vectors) != tuple(expected_vectors):
            raise ValueError(f"two-dimensional case {case_id} integer directions do not match its azimuth metadata")
        seen_directions.update(tuple(vector) for vector in vectors)
    if seen_directions != expected_directions:
        raise ValueError("two-dimensional configuration must cover all eight signed integer directions")
    if len(three_d) != 1:
        raise ValueError("the frozen Section 4.4 protocol requires exactly one three-dimensional X/Y/Z template case")
    case = three_d[0]
    case_id = str(case.get("id", ""))
    if case_id != _FORMAL_3D_CASE_ID:
        raise ValueError(f"the frozen three-dimensional template case must be '{_FORMAL_3D_CASE_ID}'")
    if [int(value) for value in case["nlag"]] != [three_dimensional_nlag, three_dimensional_nlag, three_dimensional_nlag] or [int(value) for value in case["lag_cells"]] != [1, 1, 1]:
        raise ValueError(
            f"three-dimensional case {case_id} must use frozen NLAG=({three_dimensional_nlag},{three_dimensional_nlag},{three_dimensional_nlag}) and LAG=(1,1,1)"
        )
    vectors_3d = tuple(tuple(vector) for vector in _vectors(case, 3))
    if vectors_3d != _FORMAL_3D_VECTORS:
        raise ValueError("three-dimensional production geometry must use exact X/Y/Z integer direction vectors")
    if [int(value) for value in case.get("nod", [])] != [1, 1, 1]:
        raise ValueError("three-dimensional X/Y/Z production geometry must use NOD=(1,1,1)")
    if [int(value) for value in case.get("azimuth", [])] != [0, 90, 0] or [int(value) for value in case.get("dip", [])] != [0, 0, 90]:
        raise ValueError("three-dimensional X/Y/Z production geometry must retain its frozen azimuth and dip metadata")
    if [int(value) for value in case.get("displayed_dirs", [])] != [1, 2, 3]:
        raise ValueError("three-dimensional X/Y/Z production geometry must display all three direction axes")
    expected_generators: tuple[tuple[str, str], ...] = (
        ("continuous", "continuous_iaaft_swap"),
        ("discrete", "discrete_h2_constrained_swap_walk"),
    )
    if protocol is not None:
        expected_generators = (
            _SECTION44_SURROGATE_GENERATORS[str(protocol["continuous_surrogate"])],
            _SECTION44_SURROGATE_GENERATORS[str(protocol["discrete_surrogate"])],
        )
        surrogate_spec = dict(config.get("surrogates", {}))
        expected_surrogate_fields = {
            "formal_count": int(protocol["formal_surrogates_per_dataset_property"]),
            "smoke_count": int(protocol["smoke_surrogates_per_dataset_property"]),
            "minimum_pair_difference_fraction": float(protocol["surrogate_minimum_pair_difference_fraction"]),
            "maximum_duplicate_fraction": float(protocol["surrogate_maximum_duplicate_fraction"]),
        }
        for name, expected in expected_surrogate_fields.items():
            actual = surrogate_spec.get(name)
            if actual is None or float(actual) != float(expected):
                raise ValueError(f"Section 4.4 surrogate setting {name} does not match the frozen protocol")
    for dataset in config.get("datasets", []):
        for attribute in _attribute_specs(dataset):
            generator = str(attribute["surrogate_generator"])
            field_type = str(attribute.get("field_type", ""))
            expected_field_type, expected_generator = next(
                (item for item in expected_generators if item[0] == field_type),
                ("", ""),
            )
            if not expected_field_type or generator != expected_generator:
                raise ValueError(
                    f"frozen surrogate generator for {field_type or 'unknown'} fields must be {expected_generator or 'declared'}"
                )
    confirmatory = {str(value) for value in config["confirmatory_direction_pairs"]}
    if not confirmatory or not confirmatory.issubset(case_ids):
        raise ValueError("confirmatory direction pairs must identify configured two-dimensional cases")


def _top_level_entrypoint_id(config: Mapping[str, Any], mode: str) -> str:
    """Return the run identity supplied by the only final-status owner."""

    value = config.get("top_level_entrypoint_id")
    if mode != "formal":
        return "" if value is None else str(value)
    if not isinstance(value, str) or not value.strip():
        raise ValueError("formal Section 4.4 component execution requires a top_level_entrypoint_id binding")
    return value.strip()


def run(config_path: Path, output_root: Path | None, mode: str, smoke_count_override: int | None = None) -> Path:
    config = _load(config_path)
    if mode not in {"smoke", "formal"}:
        raise ValueError("mode must be smoke or formal")
    top_level_entrypoint_id = _top_level_entrypoint_id(config, mode)
    protocol, protocol_binding = _frozen_section44_protocol(config, mode)
    sampling_mode = str(protocol["sampling_mode"]) if protocol is not None else str(config.get("template_sampling"))
    production_algorithm = str(protocol["production_algorithm"]) if protocol is not None else "AVL Based algorithm"
    audit_algorithm = str(protocol["audit_algorithm"]) if protocol is not None else "Bell algorithm"
    h2_max_relative_error = float(protocol["h2_max_relative_error"]) if protocol is not None else 0.10
    h2_weighted_relative_rmse = float(protocol["h2_weighted_relative_rmse"]) if protocol is not None else 0.05
    if sampling_mode != "exhaustive_integer_offsets":
        raise ValueError("formal Section 4.4 runner refuses non-exhaustive template sampling")
    config_dir = config_path.parent
    executable = _resolve(config_dir, str(config["template_cli"]))
    if not executable.exists():
        raise FileNotFoundError(executable)
    source_spacing = [float(value) for value in config["source_spacing"]]
    if source_spacing != [50.0, 50.0, 2.0]:
        raise ValueError("Section 4.4 formal configuration requires source_spacing=[50,50,2]")
    _validate_frozen_geometry(config, protocol)
    all_datasets = list(config["datasets"])
    datasets = all_datasets[:1] if mode == "smoke" else all_datasets
    cases_2d = list(config["directional_cases"])
    cases_3d = list(config["three_dimensional_cases"])
    if not cases_2d or not cases_3d:
        raise ValueError("at least one 2D and one 3D case is required")
    surrogate_spec = dict(config["surrogates"])
    if protocol is not None:
        if smoke_count_override is not None:
            raise ValueError("a frozen Section 4.4 protocol does not permit a surrogate-count override")
        surrogate_count = int(
            protocol["formal_surrogates_per_dataset_property"]
            if mode == "formal"
            else protocol["smoke_surrogates_per_dataset_property"]
        )
    else:
        surrogate_count = int(smoke_count_override if smoke_count_override is not None else surrogate_spec[f"{mode}_count"])
    if surrogate_count <= 0 or (mode == "formal" and surrogate_count != 199) or (mode == "smoke" and surrogate_count > 30):
        raise ValueError("formal mode requires 199 surrogates; smoke mode permits at most 30")
    root = output_root or _resolve(config_dir, str(config["output_root"]))
    run_dir = root / _timestamp()
    for name in ("configs", "inputs", "cli_outputs", "logs", "tables", "figures"):
        (run_dir / name).mkdir(parents=True, exist_ok=False)
    (run_dir / "experiment_config.json").write_bytes(config_path.read_bytes())
    if mode == "smoke":
        cases_2d, cases_3d = _smoke_cases(config, cases_2d, cases_3d)
        # The discrete H2-constrained walk must execute enough legal swaps to
        # produce a non-identical field.  This bounded count is still small
        # enough for the one-surrogate smoke run.
        surrogate_spec["swap_iterations"] = min(600, int(surrogate_spec["swap_iterations"]))
        surrogate_spec["max_attempts_per_id"] = min(12, int(surrogate_spec["max_attempts_per_id"]))
    selected_case_ids = {str(case["id"]) for case in cases_2d}
    rotation_pairs = [
        dict(pair) for pair in config["rotation_validation"]["pairs"]
        if str(pair["original_case_id"]) in selected_case_ids and str(pair["rotated_case_id"]) in selected_case_ids
    ]
    if mode == "smoke":
        expected_smoke_pairs = {
            (str(pair["original_case_id"]), str(pair["rotated_case_id"]))
            for pair in config["rotation_validation"]["pairs"]
            if bool(pair.get("include_in_smoke", False))
        }
        actual_smoke_pairs = {
            (str(pair["original_case_id"]), str(pair["rotated_case_id"]))
            for pair in rotation_pairs
        }
        if actual_smoke_pairs != expected_smoke_pairs:
            raise RuntimeError("smoke rotation matrix does not include every frozen smoke rotation pair")
    elif len(rotation_pairs) != len(config["rotation_validation"]["pairs"]):
        raise RuntimeError("formal rotation matrix does not include every frozen rotation pair")
    effective_config = {
        **config,
        "datasets": datasets,
        "directional_cases": cases_2d,
        "three_dimensional_cases": cases_3d,
        "rotation_validation": {**config["rotation_validation"], "pairs": rotation_pairs},
        "surrogates": surrogate_spec,
        "formal_section44_protocol_binding": protocol_binding,
        "top_level_entrypoint_id": top_level_entrypoint_id,
        "execution_contract": {
            "template_sampling": sampling_mode,
            "production_algorithm": production_algorithm,
            "audit_algorithm": audit_algorithm,
            "h2_max_relative_error": h2_max_relative_error,
            "h2_weighted_relative_rmse": h2_weighted_relative_rmse,
        },
        "effective_mode": mode,
        "effective_surrogate_count": surrogate_count,
    }
    (run_dir / "effective_experiment_config.json").write_text(
        json.dumps(effective_config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    native_checks = _run_native_checks(config, config_dir, run_dir / "logs")

    manifest_path, preprocessing_manifest, verified_inputs, preprocessing_contract = _validate_preprocessing_manifest(
        config, config_dir, datasets, source_spacing
    )

    confirmatory_case_ids = {str(value) for value in config["confirmatory_direction_pairs"]}
    h2_cases = [case for case in cases_2d if str(case["id"]) in confirmatory_case_ids]
    if not h2_cases:
        raise RuntimeError("no confirmatory directional templates were selected for the H2 gate")
    h2_offsets = _h2_for_cases(h2_cases)
    _write_csv(run_dir / "tables" / "H2_offsets.csv", [
        {"offset_id": index, "dy": int(offset[0]), "dx": int(offset[1]), "is_zero_lag": bool(not np.any(offset))}
        for index, offset in enumerate(h2_offsets)
    ])
    definition_rows: list[dict[str, Any]] = []
    for dimension, cases in (("2d", cases_2d), ("3d", cases_3d)):
        for case in cases:
            definition_rows.append({
                "dimension": dimension, "template_id": case["id"], "NOD": json.dumps(case["nod"]),
                "NLAG": json.dumps(case["nlag"]), "lag_cells": json.dumps(case["lag_cells"]),
                "integer_direction_vectors": json.dumps(_vectors(case, 2 if dimension == "2d" else 3)),
                "role": "confirmatory_production" if str(case["id"]) in confirmatory_case_ids else "rotation_audit",
            })
    for template_id, offsets in config.get("explicit_templates_2d", {}).items():
        definition_rows.append({"dimension": "2d", "template_id": template_id, "offsets": json.dumps(offsets), "role": "exploratory_explicit"})
    _write_csv(run_dir / "tables" / "template_definitions.csv", definition_rows)

    inventory: list[dict[str, Any]] = []
    two_d_rows: list[dict[str, Any]] = []
    three_d_rows: list[dict[str, Any]] = []
    variogram_rows: list[dict[str, Any]] = []
    acceptance_rows: list[dict[str, Any]] = []
    h2_diagnostic_rows: list[dict[str, Any]] = []
    connectivity_rows: list[dict[str, Any]] = []
    rotation_rows: list[dict[str, Any]] = []
    bell_avl_rows: list[dict[str, Any]] = []
    diversity_rows: list[dict[str, Any]] = []
    response_source: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}

    for dataset_index, dataset in enumerate(datasets):
        dataset_id = str(dataset["id"])
        for attribute_index, attribute in enumerate(_attribute_specs(dataset)):
            property_id = str(attribute["property_id"])
            values, valid, x, y, z = _read_attribute(config_dir, dataset, attribute)
            target_z = int(dataset["slice_z"])
            z_matches = np.flatnonzero(z == target_z)
            if z_matches.size != 1:
                raise RuntimeError(f"slice_z={target_z} is not present in {dataset_id}/{property_id}")
            volume_for_cli = _masked_for_cli(values, valid)
            source_input = run_dir / "inputs" / f"{dataset_id}__{property_id}__3d_original.grid"
            write_columnar_grid(source_input, volume_for_cli, str(attribute["property_name"]), x=x, y=y, z=z, allow_nan=True)
            source_path = _resolve(config_dir, str(attribute["grid_path"]))
            valid_mask_path = _resolve(config_dir, str(attribute["valid_mask_path"])) if attribute.get("valid_mask_path") else None
            inventory.append({
                "dataset_id": dataset_id, "property_id": property_id, "grid_path": str(source_path),
                "grid_sha256": sha256_file(source_path), "property_name": attribute["property_name"],
                "field_type": attribute["field_type"], "valid_mask_path": str(valid_mask_path) if valid_mask_path else "",
                "valid_mask_sha256": sha256_file(valid_mask_path) if valid_mask_path else "",
                "surrogate_generator": attribute["surrogate_generator"], "connectivity_rule": attribute["connectivity_rule"],
                "shape_zyx": json.dumps(list(values.shape)), "valid_count": int(valid.sum()),
            })
            field = values[int(z_matches[0])]
            field_valid = valid[int(z_matches[0])]
            field_for_cli = _masked_for_cli(field, field_valid)
            original_input = run_dir / "inputs" / f"{dataset_id}__{property_id}__original.grid"
            rotated_input = run_dir / "inputs" / f"{dataset_id}__{property_id}__rotated90.grid"
            write_columnar_grid(original_input, field_for_cli, str(attribute["property_name"]), x=x, y=y, allow_nan=True)
            rotated_field = np.rot90(field_for_cli, k=-1)
            rotated_valid = np.rot90(field_valid, k=-1)
            rotated_x = np.arange(1, rotated_field.shape[1] + 1, dtype=int)
            rotated_y = np.arange(1, rotated_field.shape[0] + 1, dtype=int)
            write_columnar_grid(
                rotated_input, rotated_field, str(attribute["property_name"]), x=rotated_x, y=rotated_y, allow_nan=True
            )
            global_kappa2 = _field_variance(field, field_valid)
            response_source[(dataset_id, property_id, "original", "field")] = []

            for condition, source_kind, source_input_2d, condition_field, condition_valid in (
                ("original", "original", original_input, field, field_valid),
                ("rotated90", "rotated", rotated_input, rotated_field, rotated_valid),
            ):
                for case in cases_2d:
                    target = f"{dataset_id}__{property_id}__{condition}__{case['id']}__avl"
                    result, sidecar, _ = _run_case(
                        executable=executable, run_dir=run_dir, source_grid=source_input_2d,
                        property_name=str(attribute["property_name"]), target_property=target, dataset_id=dataset_id,
                        property_id=property_id, condition=condition, case=case, dimensions=2,
                        source_spacing=source_spacing, algorithm=production_algorithm, template_sampling=sampling_mode,
                        timeout_seconds=int(config["template_cli_timeout_seconds"]),
                    )
                    records = _response_records(
                        dataset_id=dataset_id, property_id=property_id, field_type=str(attribute["field_type"]),
                        condition=condition, source_kind=source_kind, surrogate_id=None, template_id=str(case["id"]),
                        dimension="2d", algorithm=production_algorithm, result=result, rows=sidecar,
                        global_kappa2=global_kappa2, run_dir=run_dir,
                    )
                    two_d_rows.extend(records)
                    response_source[(dataset_id, property_id, condition, str(case["id"]))] = records
                    if condition == "original":
                        bell_target = f"{dataset_id}__{property_id}__original__{case['id']}__bell"
                        bell_result, bell_sidecar, _ = _run_case(
                            executable=executable, run_dir=run_dir, source_grid=source_input_2d,
                            property_name=str(attribute["property_name"]), target_property=bell_target, dataset_id=dataset_id,
                            property_id=property_id, condition="original", case=case, dimensions=2,
                            source_spacing=source_spacing, algorithm=audit_algorithm, template_sampling=sampling_mode,
                            timeout_seconds=int(config["template_cli_timeout_seconds"]),
                        )
                        bell_records = _response_records(
                            dataset_id=dataset_id, property_id=property_id, field_type=str(attribute["field_type"]),
                            condition="original", source_kind="original", surrogate_id=None, template_id=str(case["id"]),
                            dimension="2d", algorithm=audit_algorithm, result=bell_result, rows=bell_sidecar,
                            global_kappa2=global_kappa2, run_dir=run_dir,
                        )
                        two_d_rows.extend(bell_records)
                        avl_values = np.asarray([float(row["raw_c4"]) for row in records])
                        bell_values = np.asarray([float(row["raw_c4"]) for row in bell_records])
                        bell_avl_rows.append({
                            "dataset_id": dataset_id, "property_id": property_id, "dimension": "2d", "template_id": case["id"],
                            "response_count": avl_values.size,
                            "maximum_absolute_error": float(np.max(np.abs(avl_values - bell_values))),
                            "relative_l2_error": _relative_l2(avl_values, bell_values),
                        })

            generator = str(attribute["surrogate_generator"])
            surrogate_seed = int(surrogate_spec["seed"]) + dataset_index * 100_000 + attribute_index * 10_000
            if generator == "continuous_iaaft_swap":
                surrogates, audit = generate_continuous_h2_surrogates(
                    field, surrogate_count, surrogate_seed, int(surrogate_spec["iaaft_iterations"]), h2_offsets,
                    int(surrogate_spec["max_attempts_per_id"]), int(surrogate_spec["swap_iterations"]), field_valid,
                    max_absolute_relative_error=h2_max_relative_error,
                    max_weighted_relative_rmse=h2_weighted_relative_rmse,
                )
            elif generator == "discrete_h2_constrained_swap_walk":
                surrogates, audit = generate_discrete_h2_surrogates(
                    field, field_valid, surrogate_count, surrogate_seed, h2_offsets,
                    int(surrogate_spec["max_attempts_per_id"]), int(surrogate_spec["swap_iterations"]),
                    max_absolute_relative_error=h2_max_relative_error,
                    max_weighted_relative_rmse=h2_weighted_relative_rmse,
                )
            else:
                raise ValueError(f"unsupported surrogate generator {generator}")
            acceptance_rows.extend(_flatten_audit(audit, dataset_id, property_id))
            h2_diagnostic_rows.extend(_flatten_h2_diagnostics(audit, dataset_id, property_id))
            diversity = surrogate_diversity_diagnostics(surrogates)
            accepted_differences = [
                float(row["difference_fraction_from_original"])
                for row in audit
                if bool(row.get("accepted"))
            ]
            minimum_original_difference = float(min(accepted_differences)) if accepted_differences else 0.0
            diversity_rows.append({
                "dataset_id": dataset_id,
                "property_id": property_id,
                "minimum_original_difference_fraction": minimum_original_difference,
                "required_minimum_pair_difference_fraction": float(surrogate_spec["minimum_pair_difference_fraction"]),
                "maximum_duplicate_fraction": float(surrogate_spec["maximum_duplicate_fraction"]),
                **diversity,
            })
            if (
                int(diversity["unique_field_count"]) != surrogate_count
                or float(diversity["minimum_pair_difference_fraction"])
                < float(surrogate_spec["minimum_pair_difference_fraction"])
                or float(diversity["duplicate_fraction"])
                > float(surrogate_spec["maximum_duplicate_fraction"])
                or minimum_original_difference < float(surrogate_spec["minimum_original_difference_fraction"])
            ):
                raise RuntimeError(f"surrogate diversity gate failed for {dataset_id}/{property_id}")

            for surrogate_id, surrogate in enumerate(surrogates, start=1):
                condition = f"surrogate_{surrogate_id:03d}"
                surrogate_input = run_dir / "inputs" / f"{dataset_id}__{property_id}__{condition}.grid"
                write_columnar_grid(
                    surrogate_input, _masked_for_cli(surrogate, field_valid), str(attribute["property_name"]), x=x, y=y,
                    allow_nan=True,
                )
                for case in cases_2d:
                    target = f"{dataset_id}__{property_id}__{condition}__{case['id']}__avl"
                    result, sidecar, _ = _run_case(
                        executable=executable, run_dir=run_dir, source_grid=surrogate_input,
                        property_name=str(attribute["property_name"]), target_property=target, dataset_id=dataset_id,
                        property_id=property_id, condition=condition, case=case, dimensions=2,
                        source_spacing=source_spacing, algorithm=production_algorithm, template_sampling=sampling_mode,
                        timeout_seconds=int(config["template_cli_timeout_seconds"]),
                    )
                    two_d_rows.extend(_response_records(
                        dataset_id=dataset_id, property_id=property_id, field_type=str(attribute["field_type"]),
                        condition=condition, source_kind="surrogate", surrogate_id=surrogate_id, template_id=str(case["id"]),
                        dimension="2d", algorithm=production_algorithm, result=result, rows=sidecar,
                        global_kappa2=global_kappa2, run_dir=run_dir,
                    ))

            for condition, data, data_valid, kind in (("original", field, field_valid, "original"), ("rotated90", rotated_field, rotated_valid, "rotated")):
                for record in directional_semivariogram_2d(data, [1, 2, 3, 4], data_valid):
                    variogram_rows.append({"dataset_id": dataset_id, "property_id": property_id, "condition": condition, "source_kind": kind, "dimension": "2d", **record})
            for surrogate_id, surrogate in enumerate(surrogates, start=1):
                for record in directional_semivariogram_2d(surrogate, [1, 2, 3, 4], field_valid):
                    variogram_rows.append({"dataset_id": dataset_id, "property_id": property_id, "condition": f"surrogate_{surrogate_id:03d}", "source_kind": "surrogate", "surrogate_id": surrogate_id, "dimension": "2d", **record})

            for condition, data, mask, kind, sid in [("original", field, field_valid, "original", None), ("rotated90", rotated_field, rotated_valid, "rotated", None), *[(f"surrogate_{i:03d}", value, field_valid, "surrogate", i) for i, value in enumerate(surrogates, start=1)]]:
                binary, label = _mask_for_connectivity(data, mask, str(attribute["connectivity_rule"]), float(config["porosity_excursion_threshold"]))
                connectivity_rows.append({
                    "dataset_id": dataset_id, "property_id": property_id, "condition": condition, "source_kind": kind,
                    "surrogate_id": "" if sid is None else sid, "dimension": "2d", "connectivity_label": label,
                    **connectivity_metrics(binary, ("y", "x")),
                })

            # The three-dimensional path has no surrogate inference; it provides original-field description and Bell/AVL audit.
            volume_variance = _field_variance(values, valid)
            for case in cases_3d:
                target = f"{dataset_id}__{property_id}__3d__{case['id']}__avl"
                avl_result, avl_sidecar, _ = _run_case(
                    executable=executable, run_dir=run_dir, source_grid=source_input,
                    property_name=str(attribute["property_name"]), target_property=target, dataset_id=dataset_id,
                    property_id=property_id, condition="original_3d", case=case, dimensions=3,
                    source_spacing=source_spacing, algorithm=production_algorithm, template_sampling=sampling_mode,
                    timeout_seconds=int(config["template_cli_timeout_seconds"]),
                )
                avl_records = _response_records(
                    dataset_id=dataset_id, property_id=property_id, field_type=str(attribute["field_type"]), condition="original_3d",
                    source_kind="original", surrogate_id=None, template_id=str(case["id"]), dimension="3d", algorithm=production_algorithm,
                    result=avl_result, rows=avl_sidecar, global_kappa2=volume_variance,
                    run_dir=run_dir,
                )
                three_d_rows.extend(avl_records)
                bell_target = f"{dataset_id}__{property_id}__3d__{case['id']}__bell"
                bell_result, bell_sidecar, _ = _run_case(
                    executable=executable, run_dir=run_dir, source_grid=source_input,
                    property_name=str(attribute["property_name"]), target_property=bell_target, dataset_id=dataset_id,
                    property_id=property_id, condition="original_3d", case=case, dimensions=3,
                    source_spacing=source_spacing, algorithm=audit_algorithm, template_sampling=sampling_mode,
                    timeout_seconds=int(config["template_cli_timeout_seconds"]),
                )
                bell_records = _response_records(
                    dataset_id=dataset_id, property_id=property_id, field_type=str(attribute["field_type"]), condition="original_3d",
                    source_kind="original", surrogate_id=None, template_id=str(case["id"]), dimension="3d", algorithm=audit_algorithm,
                    result=bell_result, rows=bell_sidecar, global_kappa2=volume_variance,
                    run_dir=run_dir,
                )
                three_d_rows.extend(bell_records)
                avl_values = np.asarray([float(row["raw_c4"]) for row in avl_records])
                bell_values = np.asarray([float(row["raw_c4"]) for row in bell_records])
                bell_avl_rows.append({
                    "dataset_id": dataset_id, "property_id": property_id, "dimension": "3d", "template_id": case["id"],
                    "response_count": avl_values.size, "maximum_absolute_error": float(np.max(np.abs(avl_values - bell_values))),
                    "relative_l2_error": _relative_l2(avl_values, bell_values),
                })
            for record in axis_semivariogram_3d(values, [1, 2, 3], valid):
                variogram_rows.append({"dataset_id": dataset_id, "property_id": property_id, "condition": "original_3d", "source_kind": "original", "dimension": "3d", **record})
            binary_3d, label_3d = _mask_for_connectivity(values, valid, str(attribute["connectivity_rule"]), float(config["porosity_excursion_threshold"]))
            connectivity_rows.append({
                "dataset_id": dataset_id, "property_id": property_id, "condition": "original_3d", "source_kind": "original",
                "surrogate_id": "", "dimension": "3d", "connectivity_label": label_3d,
                **connectivity_metrics(binary_3d, ("z", "y", "x")),
            })

    # Each configured rotation pair must preserve the full response vector after signed-direction rotation.
    for dataset in datasets:
        dataset_id = str(dataset["id"])
        for attribute in _attribute_specs(dataset):
            property_id = str(attribute["property_id"])
            for pair in rotation_pairs:
                left = response_source[(dataset_id, property_id, "original", str(pair["original_case_id"]))]
                right = response_source[(dataset_id, property_id, "rotated90", str(pair["rotated_case_id"]))]
                left_values = np.asarray([float(row["raw_c4"]) for row in left])
                right_values = np.asarray([float(row["raw_c4"]) for row in right])
                rotation_rows.append({
                    "dataset_id": dataset_id, "property_id": property_id,
                    "original_case_id": pair["original_case_id"], "rotated_case_id": pair["rotated_case_id"],
                    "response_count": left_values.size, "shape_match": bool(left_values.shape == right_values.shape),
                    "maximum_absolute_error": float(np.max(np.abs(left_values - right_values))) if left_values.shape == right_values.shape else float("inf"),
                    "relative_l2_error": _relative_l2(left_values, right_values),
                })

    # Explicit templates are exploratory references and intentionally excluded from H2/omnibus tests.
    for dataset in datasets:
        dataset_id = str(dataset["id"])
        for attribute in _attribute_specs(dataset):
            property_id = str(attribute["property_id"])
            values, valid, x, y, z = _read_attribute(config_dir, dataset, attribute)
            index = int(np.flatnonzero(z == int(dataset["slice_z"]))[0])
            field, field_valid = values[index], valid[index]
            variance = _field_variance(field, field_valid)
            for template_id, offsets in config.get("explicit_templates_2d", {}).items():
                base = np.asarray(offsets, dtype=int)
                for angle in config["explicit_directions_degrees"]:
                    for lag in config["explicit_lags_cells"]:
                        response = template_response(np.where(field_valid, field, np.nan), transform_offsets_2d(base, int(angle), int(lag)))
                        two_d_rows.append({
                            "dataset_id": dataset_id, "property_id": property_id, "field_type": attribute["field_type"],
                            "condition": "original", "source_kind": "original", "surrogate_id": "", "template_id": template_id,
                            "direction_pair": f"angle_{angle}", "dimension": "2d", "algorithm": "reference", "calculation_path": "explicit_offset_reference",
                            "lag_cell": int(lag), "lag_positions": "", "integer_offsets": json.dumps(transform_offsets_2d(base, int(angle), int(lag)).tolist()),
                            "physical_offsets": "", "valid_translations": "", "support_count": "", "raw_c4": response,
                            "global_kappa2": variance, "standardized_k4": response / (variance * variance), "grid_x": "", "grid_y": "", "grid_z": "", "output_grid": "",
                        })

    response_summary = _response_summary(two_d_rows)
    confirmatory_cases = confirmatory_case_ids
    test_rows: list[dict[str, Any]] = []
    if surrogate_count >= 2:
        for property_id, role in (("porosity", "primary"), ("channel_occupancy", "key_secondary")):
            if not any(row["property_id"] == property_id for row in two_d_rows):
                continue
            statistic = symmetric_leave_one_out_test(_response_matrix(two_d_rows, property_id, confirmatory_cases, surrogate_count))
            test_rows.append({
                "family_id": "section44_primary" if role == "primary" else "section44_secondary",
                "hypothesis_id": f"{property_id}_omnibus_D", "role": role, "test_type": "symmetric_response_omnibus",
                "estimate": statistic["original_score"], "p_raw": statistic["empirical_p_value"],
                "active_cell_count": statistic["active_cell_count"], "constant_cell_count": statistic["constant_cell_count"],
                "numerical_floor": statistic["numerical_floor"], "inference_status": "evaluated",
            })
        if any(row["property_id"] == "channel_occupancy" and row["dimension"] == "2d" for row in connectivity_rows):
            statistic = symmetric_leave_one_out_test(_connectivity_test_matrix(connectivity_rows, "channel_occupancy", surrogate_count))
            test_rows.append({
                "family_id": "section44_secondary", "hypothesis_id": "channel_strict_majority_largest_component_D",
                "role": "key_secondary", "test_type": "symmetric_connectivity_omnibus", "estimate": statistic["original_score"],
                "p_raw": statistic["empirical_p_value"], "active_cell_count": statistic["active_cell_count"],
                "constant_cell_count": statistic["constant_cell_count"], "numerical_floor": statistic["numerical_floor"],
                "inference_status": "evaluated",
            })
        secondary = [row for row in test_rows if row["family_id"] == "section44_secondary"]
        adjusted = {str(row["hypothesis_id"]): row for row in holm_step_down(secondary)}
        for row in test_rows:
            if row["family_id"] == "section44_secondary":
                row.update({key: value for key, value in adjusted[str(row["hypothesis_id"])].items() if key in {"holm_rank", "p_holm", "reject_holm_0_05"}})
            else:
                row.update({"holm_rank": "", "p_holm": "", "reject_holm_0_05": ""})
    else:
        # A one-surrogate smoke run validates the computational chain but has
        # fewer than three exchangeable candidates, so it cannot supply an
        # omnibus surrogate p-value or a Holm adjustment.
        test_rows.append({
            "family_id": "smoke_only",
            "hypothesis_id": "omnibus_not_evaluated",
            "role": "not_for_inference",
            "test_type": "not_evaluated_single_surrogate_smoke",
            "estimate": "",
            "p_raw": "",
            "active_cell_count": "",
            "constant_cell_count": "",
            "numerical_floor": "",
            "holm_rank": "",
            "p_holm": "",
            "reject_holm_0_05": "",
            "inference_status": "not_evaluated_single_surrogate_smoke",
        })

    _write_csv(run_dir / "tables" / "dataset_inventory.csv", inventory)
    _write_csv(run_dir / "tables" / "template_2d_response_long.csv", two_d_rows)
    _write_csv(run_dir / "tables" / "template_3d_response_long.csv", three_d_rows)
    _write_csv(run_dir / "tables" / "directional_variogram_long.csv", variogram_rows)
    _write_csv(run_dir / "tables" / "surrogate_acceptance.csv", acceptance_rows)
    _write_csv(run_dir / "tables" / "surrogate_h2_diagnostics.csv", h2_diagnostic_rows)
    _write_csv(run_dir / "tables" / "surrogate_diversity.csv", diversity_rows)
    _write_csv(run_dir / "tables" / "surrogate_response_summary.csv", response_summary)
    _write_csv(run_dir / "tables" / "rotation_validation_long.csv", rotation_rows)
    _write_csv(run_dir / "tables" / "bell_avl_accuracy.csv", bell_avl_rows)
    _write_csv(run_dir / "tables" / "connectivity_metrics.csv", connectivity_rows)
    _write_csv(run_dir / "tables" / "multiplicity_decisions.csv", test_rows)
    figures = _render_figures_from_tables(run_dir)

    tolerance = {
        "max_absolute_error": float(protocol["bell_avl_max_abs"])
        if protocol is not None else float(config["algorithm_audit"]["max_absolute_error"]),
        "max_relative_l2_error": float(protocol["bell_avl_relative_l2"])
        if protocol is not None else float(config["algorithm_audit"]["max_relative_l2_error"]),
    }
    rotation_relative_l2 = (
        float(protocol["rotation_relative_l2"])
        if protocol is not None else float(config["rotation_validation"]["max_relative_l2_error"])
    )
    expected_surrogate_keys = {
        (str(dataset["id"]), str(attribute["property_id"]), surrogate_id)
        for dataset in datasets
        for attribute in _attribute_specs(dataset)
        for surrogate_id in range(1, surrogate_count + 1)
    }
    accepted_surrogate_rows = [row for row in acceptance_rows if bool(row.get("accepted"))]
    gates = {
        "template_sampling_is_exhaustive_integer_offsets": True,
        "frozen_section44_protocol_binding": protocol_binding is not None or mode == "smoke",
        "preprocessing_manifest_inputs_verified": bool(verified_inputs),
        "preprocessing_axis_mapping_and_fixture_verified": bool(preprocessing_contract),
        "native_core_checks": all(check["returncode"] == 0 for check in native_checks),
        "surrogate_count_complete": _accepted_surrogate_keys(acceptance_rows) == expected_surrogate_keys,
        # Rejected candidates stay in the audit table. Only accepted fields
        # enter template response calculations and this acceptance gate.
        "surrogate_h2_gate": len(accepted_surrogate_rows) == len(expected_surrogate_keys),
        "surrogate_diversity": all(
            float(row["duplicate_fraction"]) <= float(row["maximum_duplicate_fraction"])
            and float(row["minimum_pair_difference_fraction"])
            >= float(row["required_minimum_pair_difference_fraction"])
            and float(row["minimum_original_difference_fraction"])
            >= float(surrogate_spec["minimum_original_difference_fraction"])
            for row in diversity_rows
        ),
        "bell_avl_max_absolute_error": all(float(row["maximum_absolute_error"]) <= float(tolerance["max_absolute_error"]) for row in bell_avl_rows),
        "bell_avl_relative_l2_error": all(float(row["relative_l2_error"]) <= float(tolerance["max_relative_l2_error"]) for row in bell_avl_rows),
        "rotation_signed_direction_max_absolute_error": all(
            float(row["maximum_absolute_error"]) <= float(config["rotation_validation"]["max_absolute_error"])
            for row in rotation_rows
        ),
        "rotation_signed_direction_relative_l2_error": all(
            float(row["relative_l2_error"]) <= rotation_relative_l2
            for row in rotation_rows
        ),
        "required_figures": len(figures) >= 3,
    }
    # This component cannot certify the multi-section study on its own.  The
    # top-level formal entrypoint is the only owner of a final PASS state.
    status = ("SMOKE_PASS" if mode == "smoke" else "DIAGNOSTIC") if all(gates.values()) else "FAIL"
    gate = {
        "status": status,
        "mode": mode,
        "top_level_entrypoint_id": top_level_entrypoint_id,
        "gates": gates,
        "surrogate_count": surrogate_count,
        "formal_section44_protocol_sha256": protocol_binding["sha256"] if protocol_binding is not None else "",
    }
    (run_dir / "gate.json").write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = _write_report(run_dir, status, gates, surrogate_count)
    index = _write_component_index(run_dir)
    hash_rows = _hash_table(run_dir)
    _write_csv(run_dir / "sha256_manifest.csv", hash_rows)
    manifest = {
        "experiment": "template_directional_20261010", "status": status, "mode": mode,
        "top_level_entrypoint_id": top_level_entrypoint_id,
        "created_utc": datetime.now(timezone.utc).isoformat(), "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "effective_config": {
            "path": str(run_dir / "effective_experiment_config.json"),
            "sha256": sha256_file(run_dir / "effective_experiment_config.json"),
        },
        "formal_section44_protocol_binding": protocol_binding,
        "execution_contract": effective_config["execution_contract"],
        "runner": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__).resolve())},
        "template_cli": {"path": str(executable), "sha256": sha256_file(executable)},
        "preprocessing_manifest": {
            "path": str(manifest_path), "sha256": sha256_file(manifest_path),
            "schema_version": preprocessing_manifest["schema_version"],
        },
        "preprocessing_contract": preprocessing_contract,
        "verified_preprocessing_inputs": verified_inputs,
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
        "gates": gates, "reports": [report.name, index.name],
        "native_checks": [{"name": row["name"], "command": row["command"], "stdout": str(row["stdout"]), "stderr": str(row["stderr"])} for row in native_checks],
        "outputs": hash_rows,
    }
    (run_dir / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return run_dir


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "experiment_config.json")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--mode", choices=("smoke", "formal"), default="smoke")
    parser.add_argument("--smoke", action="store_true", help="compatibility alias for --mode smoke")
    parser.add_argument("--smoke-count", type=int, help="development override capped at 30 accepted surrogate fields")
    args = parser.parse_args(argv)
    mode = "smoke" if args.smoke else args.mode
    print(run(args.config.resolve(), args.output_root.resolve() if args.output_root else None, mode, args.smoke_count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
