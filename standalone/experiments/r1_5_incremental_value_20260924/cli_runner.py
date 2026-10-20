from __future__ import annotations

import hashlib
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from grid_io import GridData, align_maps, read_cli_output


ATOL = 1.0e-12
RTOL = 1.0e-10
CUMULANT_CENTERING_VERSION = "order_ge_2.v1"
FORMAL_WEIGHT_MODES = frozenset({"uniform", "coordinate_gaussian", "soft_high_porosity"})
WINDOW_SUPPORT_COLUMNS = (
    "x",
    "y",
    "z",
    "weight_sum",
    "weight_squared_sum",
    "effective_sample_size",
    "sample_count",
    "map_included",
    "valid",
    "valid_reason",
)


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bool_value(value: object, fallback: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return fallback
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "on"}:
            return True
        if normalized in {"0", "false", "no", "off"}:
            return False
        raise ValueError(f"Invalid Boolean value: {value!r}")
    if isinstance(value, (int, np.integer)):
        return bool(value)
    raise ValueError(f"Invalid Boolean value: {value!r}")


def _safe_token(value: object, *, field: str) -> str:
    token = str(value).strip()
    if not token:
        raise ValueError(f"{field} must not be empty")
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "-", token).strip(".-")
    if not normalized:
        raise ValueError(f"{field} does not contain a usable path-safe token: {value!r}")
    return normalized


def _window_label(window_size: tuple[int, int, int]) -> str:
    return "x".join(str(int(value)) for value in window_size)


def _read_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            metadata[key.strip()] = value.strip()
    return metadata


def _parse_sidecar_bool(series: pd.Series, column: str) -> pd.Series:
    normalized = series.astype(str).str.strip().str.lower()
    values = normalized.map({"1": True, "true": True, "yes": True, "on": True, "0": False, "false": False, "no": False, "off": False})
    if values.isna().any():
        invalid = series[values.isna()].head().tolist()
        raise ValueError(f"Invalid {column} value(s) in window support sidecar: {invalid}")
    return values.astype(bool)


def read_window_support_sidecar(path: str | Path) -> pd.DataFrame:
    """Read and validate the C++ per-window support sidecar.

    The sidecar is intentionally retained at attempted-window granularity.  A
    production map only contains the rows marked ``map_included``; callers use
    ``valid`` and ``valid_reason`` to define the formal common support.
    """

    sidecar_path = Path(path)
    if not sidecar_path.is_file():
        raise FileNotFoundError(f"C++ window support sidecar does not exist: {sidecar_path}")
    frame = pd.read_csv(sidecar_path)
    missing = [column for column in WINDOW_SUPPORT_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Window support sidecar is missing required columns {missing}: {sidecar_path}")
    result = frame.loc[:, list(WINDOW_SUPPORT_COLUMNS)].copy()
    for column in ("x", "y", "z", "sample_count"):
        result[column] = pd.to_numeric(result[column], errors="raise").astype(int)
    for column in ("weight_sum", "weight_squared_sum", "effective_sample_size"):
        result[column] = pd.to_numeric(result[column], errors="raise").astype(float)
    result["map_included"] = _parse_sidecar_bool(result["map_included"], "map_included")
    result["valid"] = _parse_sidecar_bool(result["valid"], "valid")
    result["valid_reason"] = result["valid_reason"].fillna("").astype(str)
    if (result["valid_reason"].str.len() == 0).any():
        raise ValueError(f"Window support sidecar contains an empty valid_reason: {sidecar_path}")
    if result.duplicated(["x", "y", "z"]).any():
        duplicate = result.loc[result.duplicated(["x", "y", "z"], keep=False), ["x", "y", "z"]].head()
        raise ValueError(f"Window support sidecar has duplicate centers: {duplicate.to_dict('records')}")
    return result.sort_values(["z", "y", "x"]).reset_index(drop=True)


def validate_cli_window_alignment(
    map_frame: pd.DataFrame,
    support_frame: pd.DataFrame,
    expected_coordinates: pd.DataFrame | None = None,
    *,
    require_valid_support: bool = False,
) -> dict[str, int | str | bool]:
    """Verify map, sidecar, and optionally ``WindowBatch`` centers by coordinates."""

    coordinate_columns = ["x", "y", "z"]
    required_map = set(coordinate_columns + ["value"])
    if not required_map.issubset(map_frame.columns):
        raise ValueError("C++ map must contain x, y, z, and value columns")
    if not set(WINDOW_SUPPORT_COLUMNS).issubset(support_frame.columns):
        raise ValueError("C++ support frame does not contain the required sidecar schema")
    if map_frame.duplicated(coordinate_columns).any():
        raise ValueError("C++ map has duplicate output coordinates")
    included = support_frame.loc[support_frame["map_included"], coordinate_columns]
    if included.duplicated(coordinate_columns).any():
        raise ValueError("C++ support sidecar has duplicate included coordinates")
    map_coordinates = map_frame.loc[:, coordinate_columns]
    map_sidecar = map_coordinates.merge(included, on=coordinate_columns, how="outer", indicator=True)
    if not (map_sidecar["_merge"] == "both").all():
        counts = map_sidecar["_merge"].value_counts().to_dict()
        raise ValueError(f"C++ map and support-sidecar coordinates differ: {counts}")
    if require_valid_support:
        invalid = support_frame.loc[support_frame["map_included"] & ~support_frame["valid"]]
        if not invalid.empty:
            reasons = invalid["valid_reason"].value_counts().to_dict()
            raise ValueError(f"C++ map includes invalid support rows: {reasons}")
    expected_count = None
    if expected_coordinates is not None:
        if not set(coordinate_columns).issubset(expected_coordinates.columns):
            raise ValueError("Expected WindowBatch coordinates must contain x, y, z")
        expected = expected_coordinates.loc[:, coordinate_columns].copy()
        for column in coordinate_columns:
            expected[column] = pd.to_numeric(expected[column], errors="raise").astype(int)
        if expected.duplicated(coordinate_columns).any():
            raise ValueError("Expected WindowBatch coordinates are not unique")
        expected_count = int(len(expected))
        comparison = map_coordinates.merge(expected, on=coordinate_columns, how="outer", indicator=True)
        if not (comparison["_merge"] == "both").all():
            counts = comparison["_merge"].value_counts().to_dict()
            raise ValueError(f"C++ map and WindowBatch center coordinates differ: {counts}")
    return {
        "passed": True,
        "status": "PASS",
        "map_node_count": int(len(map_frame)),
        "sidecar_row_count": int(len(support_frame)),
        "included_node_count": int(len(included)),
        "expected_node_count": expected_count if expected_count is not None else int(len(map_frame)),
    }


@dataclass(frozen=True)
class CliResult:
    route: str
    config_path: Path
    output_path: Path
    metadata_path: Path
    stdout_path: Path
    stderr_path: Path
    frame: pd.DataFrame
    case_id: str = ""
    target_property: str = ""
    sidecar_path: Path | None = None
    sidecar_frame: pd.DataFrame | None = None
    authoritative_cplusplus: bool = True
    elapsed_seconds: float = float("nan")


@dataclass(frozen=True)
class CppFeatureMapSpec:
    """One formal C++ cumulant-map request.

    ``grid`` is the immutable source grid for the requested dataset.  The
    public matrix API accepts only the three formal weighting schemes and
    assigns AVL to uniform maps and the weighted-moment route to nonuniform
    maps.  ``parameter_role`` distinguishes frozen primary settings from the
    prespecified sensitivity settings without changing the map identity.
    """

    dataset_id: str
    grid: GridData
    window_size: tuple[int, int, int]
    order: int
    weight_config: Mapping[str, Any]
    section: str | None = None
    weighting_id: str | None = None
    parameter_set: str | None = None
    parameter_role: str = "primary"
    move_step: tuple[int, int, int] = (1, 1, 1)
    source_spacing: tuple[float, float, float] = (1.0, 1.0, 1.0)


@dataclass(frozen=True)
class CppFeatureMapArtifact:
    """One coordinate-validated C++ map and its per-window support sidecar."""

    spec: CppFeatureMapSpec
    cli_result: CliResult
    map_frame: pd.DataFrame
    support_frame: pd.DataFrame
    algorithm: str
    coordinate_validation: str
    manifest_row: Mapping[str, Any]
    authoritative_cplusplus: bool = True


@dataclass(frozen=True)
class CppFeatureMatrixResult:
    """Artifacts and frozen manifest returned by :func:`run_cpp_feature_matrix`."""

    artifacts: tuple[CppFeatureMapArtifact, ...]
    manifest: pd.DataFrame
    manifest_path: Path
    uniform_manifest_path: Path
    weighted_manifest_path: Path


def _as_csv(values) -> str:
    return ",".join(str(value) for value in values)


def _experiment_weight_switches(weight_config: Mapping[str, Any]) -> tuple[str, dict[str, bool]]:
    """Resolve the three schemes admitted to the R1-5 experiment chain."""
    aliases = {
        "coordinate": "coordinate_gaussian",
        "spatial_gaussian": "coordinate_gaussian",
        "high_porosity": "soft_high_porosity",
    }
    mode = aliases.get(str(weight_config.get("mode", "uniform")).lower(), str(weight_config.get("mode", "uniform")).lower())
    if mode not in FORMAL_WEIGHT_MODES:
        raise ValueError(f"R1-5 experiments support uniform, coordinate_gaussian, and soft_high_porosity; got {mode}")
    switches = {
        "enable_coordinate_gaussian": mode == "coordinate_gaussian",
        "enable_soft_high_porosity": mode == "soft_high_porosity",
    }
    for key, active in switches.items():
        if key in weight_config and _bool_value(weight_config[key]) != active:
            raise ValueError(f"{key} conflicts with weighting mode {mode}")
    if (
        _bool_value(weight_config.get("enable_fixed_value_gaussian", False))
        or _bool_value(weight_config.get("enable_adaptive_upper_tail", False))
        or _bool_value(weight_config.get("enable_joint_geometric_mean", False))
    ):
        raise ValueError(
            "fixed_value_gaussian, adaptive_upper_tail, and joint_geometric_mean "
            "are excluded from the R1-5 experiment chain"
        )
    if not _bool_value(weight_config.get("normalize", True), True):
        raise ValueError("R1-5 cumulant experiments require normalized weights")
    return mode, switches


def build_ini(
    source_grid: Path,
    source_property: str,
    target_property: str,
    output_dir: Path,
    window_size: tuple[int, int, int],
    move_step: tuple[int, int, int],
    source_origin: tuple[int, int, int],
    source_spacing: tuple[float, float, float],
    coordinate_dimensions: int,
    algorithm: str,
    weighted: bool,
    weight_config: Mapping[str, Any],
    order: int = 4,
) -> str:
    mode, switches = _experiment_weight_switches(weight_config)
    coordinate_center = ""
    if mode == "coordinate_gaussian":
        coordinate_center = "weight_center_mode=window_center\n"
    zero_first_order_moment = _bool_value(weight_config.get("zero_first_order_moment", False))
    return f"""[general]
calculate_target=Cumulant
calculate_algorithm={algorithm}
source_grid={source_grid.as_posix()}
source_property={source_property}
coordinate_dimensions={coordinate_dimensions}
target_property={target_property}
source_origin={_as_csv(source_origin)}
source_spacing={_as_csv(source_spacing)}

[sliding_window]
window_size={_as_csv(window_size)}
move_step={_as_csv(move_step)}
moment_cumulant_order={order}
center_cumulant_input_for_order_ge_2=true
require_full_window=true
require_complete_window=true
require_anchor_present=true
skip_missing_samples=false
output_at_window_center=true
invert_window_output=false
enable_percentile_clipping=false
minimum_samples=1

[weighted_moment]
enable_weighted_moment={'true' if weighted else 'false'}
weight_mode={mode}
enable_coordinate_gaussian={'true' if switches['enable_coordinate_gaussian'] else 'false'}
enable_soft_high_porosity={'true' if switches['enable_soft_high_porosity'] else 'false'}
bandwidth={weight_config.get('bandwidth', 1.0)}
min_weight={weight_config.get('min_weight', 0.0)}
normalize_weights={'true' if _bool_value(weight_config.get('normalize', True), True) else 'false'}
zero_first_order_moment={'true' if zero_first_order_moment else 'false'}
high_porosity_threshold={weight_config.get('high_porosity_threshold', 0.25)}
soft_transition_width={weight_config.get('soft_transition_width', 0.02)}
soft_weight_strength={weight_config.get('soft_weight_strength', 1.0)}
{coordinate_center}

[runtime]
output_dir={output_dir.as_posix()}
write_window_support_sidecar=true
"""


def run_cli_case(
    executable: Path,
    grid: GridData,
    dataset_id: str,
    window_size: tuple[int, int, int],
    move_step: tuple[int, int, int],
    route: str,
    weight_config: Mapping[str, Any],
    source_spacing: tuple[float, float, float],
    run_dir: Path,
    timeout_seconds: float,
    order: int = 4,
    *,
    case_id: str | None = None,
    target_property: str | None = None,
    require_window_support_sidecar: bool = True,
) -> CliResult:
    """Run one C++ moving-window map and load its support sidecar.

    Existing callers can continue using the positional interface.  Formal map
    production supplies a stable ``case_id`` so every INI, log, map, metadata
    file, and sidecar has a collision-free provenance key.
    """

    executable = Path(executable).resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"movewindow_cli executable does not exist: {executable}")
    if int(order) < 1:
        raise ValueError("Cumulant order must be positive")
    weighted = route == "weighted"
    if route not in {"bell", "avl", "weighted"}:
        raise ValueError(f"Unsupported C++ route: {route}")
    algorithm = "Bell algorithm" if route == "bell" else "AVL Based algorithm"
    if weighted:
        algorithm = "Weighted moment"
    window_label = "w" + "x".join(str(value) for value in window_size)
    resolved_case_id = _safe_token(case_id or f"{dataset_id}_{window_label}_r{order}_{route}", field="case_id")
    resolved_target_property = _safe_token(target_property or f"cumulant{order}_{route}", field="target_property")
    run_dir = Path(run_dir)
    config_path = run_dir / "configs" / f"{resolved_case_id}.ini"
    output_dir = run_dir / "maps" / "cpp" / resolved_case_id
    logs_dir = run_dir / "logs"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    origin = (int(np.min(grid.x)), int(np.min(grid.y)), int(np.min(grid.z)))
    coordinate_dimensions = 2 if grid.values.shape[0] == 1 else 3
    config_path.write_text(
        build_ini(
            grid.source_path,
            grid.property_name,
            resolved_target_property,
            output_dir,
            window_size,
            move_step,
            origin,
            source_spacing,
            coordinate_dimensions,
            algorithm,
            weighted,
            weight_config,
            int(order),
        ),
        encoding="utf-8",
    )
    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable), str(config_path)],
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    elapsed_seconds = time.perf_counter() - started
    stdout_path = logs_dir / f"{resolved_case_id}.stdout.txt"
    stderr_path = logs_dir / f"{resolved_case_id}.stderr.txt"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"CLI failed for {resolved_case_id} (code {completed.returncode}); see {stderr_path}")
    output_path = output_dir / f"{resolved_target_property}.txt"
    metadata_path = output_dir / f"{resolved_target_property}.meta.txt"
    if not output_path.exists() or not metadata_path.exists():
        raise RuntimeError(f"CLI did not produce the expected files for {resolved_case_id}")
    sidecar_path = output_dir / f"{resolved_target_property}.window_support.csv"
    if require_window_support_sidecar and not sidecar_path.is_file():
        raise RuntimeError(f"CLI did not produce the required window support sidecar for {resolved_case_id}")
    frame = read_cli_output(output_path)
    sidecar_frame = read_window_support_sidecar(sidecar_path) if sidecar_path.is_file() else None
    if sidecar_frame is not None:
        validate_cli_window_alignment(frame, sidecar_frame)
    return CliResult(
        route=route,
        config_path=config_path,
        output_path=output_path,
        metadata_path=metadata_path,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        frame=frame,
        case_id=resolved_case_id,
        target_property=resolved_target_property,
        sidecar_path=sidecar_path if sidecar_path.is_file() else None,
        sidecar_frame=sidecar_frame,
        authoritative_cplusplus=True,
        elapsed_seconds=float(elapsed_seconds),
    )


def _default_parameter_set(mode: str, weight_config: Mapping[str, Any]) -> str:
    if mode == "uniform":
        return "uniform"
    if mode == "coordinate_gaussian":
        return f"coordinate_gaussian_sigma_{float(weight_config.get('bandwidth', 1.0)):.2f}"
    if mode == "soft_high_porosity":
        return (
            "soft_high_porosity_"
            f"threshold_{float(weight_config.get('high_porosity_threshold', 0.25)):.2f}_"
            f"width_{float(weight_config.get('soft_transition_width', 0.02)):.2f}_"
            f"strength_{float(weight_config.get('soft_weight_strength', 1.0)):.2f}"
        )
    raise ValueError(f"Unsupported formal weighting mode: {mode}")


def _resolve_feature_spec(spec: CppFeatureMapSpec) -> dict[str, Any]:
    weight_config = dict(spec.weight_config)
    mode, _ = _experiment_weight_switches(weight_config)
    if not 1 <= int(spec.order) <= 6:
        raise ValueError(f"Formal C++ cumulant-map order must lie in [1, 6], got {spec.order}")
    window_size = tuple(int(value) for value in spec.window_size)
    move_step = tuple(int(value) for value in spec.move_step)
    source_spacing = tuple(float(value) for value in spec.source_spacing)
    if len(window_size) != 3 or any(value <= 0 for value in window_size):
        raise ValueError(f"window_size must contain three positive values, got {spec.window_size}")
    if len(move_step) != 3:
        raise ValueError(f"move_step must contain three values, got {spec.move_step}")
    if len(source_spacing) != 3 or any(value <= 0.0 for value in source_spacing):
        raise ValueError("source_spacing must contain three positive values")
    if spec.grid.values.shape[0] == 1:
        if move_step[0] <= 0 or move_step[1] <= 0 or move_step[2] != 0:
            raise ValueError("2-D formal maps require positive X/Y moves and zero Z move")
    elif any(value <= 0 for value in move_step):
        raise ValueError("3-D formal maps require three positive move steps")
    role = str(spec.parameter_role).strip().lower()
    if role not in {"primary", "sensitivity"}:
        raise ValueError(f"parameter_role must be 'primary' or 'sensitivity', got {spec.parameter_role!r}")
    if role == "sensitivity" and mode == "uniform":
        raise ValueError("Uniform belongs to the primary C++ map matrix; sensitivity maps must be nonuniform")
    section = "4.5" if mode == "uniform" else "4.6"
    if spec.section is not None and str(spec.section).strip() not in {"", section}:
        raise ValueError(f"{mode} C++ maps belong to Section {section}, not {spec.section!r}")
    weighting_id = _safe_token(spec.weighting_id or mode, field="weighting_id")
    parameter_set = _safe_token(spec.parameter_set or _default_parameter_set(mode, weight_config), field="parameter_set")
    dataset_id = _safe_token(spec.dataset_id, field="dataset_id")
    algorithm = "AVL" if mode == "uniform" else "Weighted moment"
    route = "avl" if mode == "uniform" else "weighted"
    return {
        "dataset_id": dataset_id,
        "mode": mode,
        "weight_config": weight_config,
        "window_size": window_size,
        "move_step": move_step,
        "source_spacing": source_spacing,
        "section": section,
        "weighting_id": weighting_id,
        "parameter_set": parameter_set,
        "parameter_role": role,
        "algorithm": algorithm,
        "route": route,
    }


def _feature_case_id(specification: Mapping[str, Any], order: int) -> str:
    readable = "_".join(
        (
            "cpp",
            str(specification["dataset_id"]),
            "w" + _window_label(specification["window_size"]),
            str(specification["weighting_id"]),
            str(specification["parameter_set"]),
            f"kappa{int(order)}",
        )
    )
    digest = hashlib.sha256(readable.encode("utf-8")).hexdigest()[:12]
    return _safe_token(f"{readable[:112]}_{digest}", field="case_id")


def build_cpp_feature_map_specs(
    datasets: Mapping[str, GridData],
    windows: Sequence[tuple[int, int, int]],
    primary_weights: Sequence[Mapping[str, Any]],
    sensitivity_weights: Sequence[Mapping[str, Any]] = (),
    *,
    orders: Sequence[int] = (1, 2, 3, 4, 5, 6),
    move_step: tuple[int, int, int] = (1, 1, 1),
    source_spacing: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> tuple[CppFeatureMapSpec, ...]:
    """Build the deterministic formal-map cartesian product without running it.

    The caller supplies the already-loaded source grids.  Primary rows retain
    the three main weighting IDs; sensitivity rows retain their explicit
    configuration IDs.  The resulting tuple can be passed directly to
    :func:`run_cpp_feature_matrix` for a full run or a deliberately smaller
    smoke subset.
    """

    if not datasets:
        raise ValueError("At least one dataset is required to build C++ map specifications")
    if not windows:
        raise ValueError("At least one window is required to build C++ map specifications")
    if not primary_weights:
        raise ValueError("At least one primary weighting configuration is required")
    resolved_orders = tuple(int(order) for order in orders)
    if not resolved_orders or any(order < 1 or order > 6 for order in resolved_orders):
        raise ValueError("Formal C++ map orders must be a nonempty subset of [1, 6]")

    specs: list[CppFeatureMapSpec] = []
    for parameter_role, weightings in (("primary", primary_weights), ("sensitivity", sensitivity_weights)):
        for weight_mapping in weightings:
            mapping = dict(weight_mapping)
            mode, _ = _experiment_weight_switches(mapping)
            parameter_set = _safe_token(mapping.get("id", _default_parameter_set(mode, mapping)), field="weighting id")
            for dataset_id, grid in datasets.items():
                for window_size in windows:
                    for order in resolved_orders:
                        specs.append(
                            CppFeatureMapSpec(
                                dataset_id=str(dataset_id),
                                grid=grid,
                                window_size=tuple(int(value) for value in window_size),
                                order=order,
                                weight_config=mapping,
                                section="4.5" if mode == "uniform" else "4.6",
                                weighting_id=mode,
                                parameter_set=parameter_set,
                                parameter_role=parameter_role,
                                move_step=move_step,
                                source_spacing=source_spacing,
                            )
                        )
    return tuple(specs)


def build_direct_cpp_feature_map_specs(
    datasets: Mapping[str, GridData],
    *,
    section45_windows: Sequence[tuple[int, int, int]],
    section46_primary_window: tuple[int, int, int],
    uniform_weight: Mapping[str, Any],
    nonuniform_primary_weights: Sequence[Mapping[str, Any]],
    sensitivity_weights: Sequence[Mapping[str, Any]],
    orders: Sequence[int] = (1, 2, 3, 4, 5, 6),
    sensitivity_orders: Sequence[int] = (2, 6),
    sensitivity_dataset_id: str = "z45_56",
    move_step: tuple[int, int, int] = (1, 1, 1),
    source_spacing: tuple[float, float, float] = (50.0, 50.0, 2.0),
) -> tuple[CppFeatureMapSpec, ...]:
    """Build the frozen 90-primary plus 12-sensitivity direct-statistics matrix.

    The direct protocol is intentionally asymmetric.  Section 4.5 evaluates
    uniform windows at three scales, Section 4.6 holds the primary scale fixed
    while changing the two nonuniform weights, and the six sensitivity settings
    are restricted to the representative subvolume and orders 2 and 6.
    """

    if not datasets:
        raise ValueError("At least one dataset is required for the direct map protocol")
    resolved_orders = tuple(int(order) for order in orders)
    resolved_sensitivity_orders = tuple(int(order) for order in sensitivity_orders)
    if set(resolved_orders) != {1, 2, 3, 4, 5, 6}:
        raise ValueError("Direct primary maps require exactly orders 1 through 6")
    if set(resolved_sensitivity_orders) != {2, 6}:
        raise ValueError("Direct sensitivity maps require exactly orders 2 and 6")
    if sensitivity_dataset_id not in datasets:
        raise ValueError(f"Direct sensitivity dataset is unavailable: {sensitivity_dataset_id}")
    if len(section45_windows) != 3:
        raise ValueError("Direct Section 4.5 requires exactly three window sizes")

    uniform_mode, _ = _experiment_weight_switches(uniform_weight)
    if uniform_mode != "uniform":
        raise ValueError("Direct Section 4.5 requires a uniform weighting configuration")
    nonuniform = [dict(weight) for weight in nonuniform_primary_weights]
    nonuniform_modes = {_experiment_weight_switches(weight)[0] for weight in nonuniform}
    if nonuniform_modes != {"coordinate_gaussian", "soft_high_porosity"}:
        raise ValueError("Direct Section 4.6 requires coordinate_gaussian and soft_high_porosity exactly once")
    if len(nonuniform) != 2:
        raise ValueError("Direct Section 4.6 requires exactly two nonuniform primary weight settings")
    if len(sensitivity_weights) != 6:
        raise ValueError("Direct protocol requires exactly six prespecified sensitivity settings")

    specs: list[CppFeatureMapSpec] = []
    for dataset_id, grid in datasets.items():
        for window_size in section45_windows:
            for order in resolved_orders:
                specs.append(
                    CppFeatureMapSpec(
                        dataset_id=str(dataset_id),
                        grid=grid,
                        window_size=tuple(int(value) for value in window_size),
                        order=order,
                        weight_config=dict(uniform_weight),
                        section="4.5",
                        weighting_id="uniform",
                        parameter_set=str(uniform_weight.get("id", "uniform")),
                        parameter_role="primary",
                        move_step=move_step,
                        source_spacing=source_spacing,
                    )
                )
        for weight in nonuniform:
            mode, _ = _experiment_weight_switches(weight)
            parameter_set = str(weight.get("id", _default_parameter_set(mode, weight)))
            for order in resolved_orders:
                specs.append(
                    CppFeatureMapSpec(
                        dataset_id=str(dataset_id),
                        grid=grid,
                        window_size=tuple(int(value) for value in section46_primary_window),
                        order=order,
                        weight_config=weight,
                        section="4.6",
                        weighting_id=mode,
                        parameter_set=parameter_set,
                        parameter_role="primary",
                        move_step=move_step,
                        source_spacing=source_spacing,
                    )
                )
    sensitivity_grid = datasets[sensitivity_dataset_id]
    for weight in sensitivity_weights:
        mapping = dict(weight)
        mode, _ = _experiment_weight_switches(mapping)
        if mode == "uniform":
            raise ValueError("Direct sensitivity maps must use a nonuniform weighting setting")
        parameter_set = str(mapping.get("id", _default_parameter_set(mode, mapping)))
        for order in resolved_sensitivity_orders:
            specs.append(
                CppFeatureMapSpec(
                    dataset_id=sensitivity_dataset_id,
                    grid=sensitivity_grid,
                    window_size=tuple(int(value) for value in section46_primary_window),
                    order=order,
                    weight_config=mapping,
                    section="4.6",
                    weighting_id=mode,
                    parameter_set=parameter_set,
                    parameter_role="sensitivity",
                    move_step=move_step,
                    source_spacing=source_spacing,
                )
            )
    return tuple(specs)


def validate_direct_cpp_feature_map_protocol(
    specs: Sequence[CppFeatureMapSpec],
    *,
    expected_primary_count: int = 90,
    expected_sensitivity_count: int = 12,
) -> dict[str, int]:
    """Validate every direct-protocol identity before a C++ production run."""

    resolved = [dict(_resolve_feature_spec(spec), order=int(spec.order)) for spec in specs]
    keys = [
        (
            item["section"],
            item["dataset_id"],
            _window_label(item["window_size"]),
            item["weighting_id"],
            item["parameter_set"],
            item["parameter_role"],
            item["order"],
            item["algorithm"],
        )
        for item in resolved
    ]
    if len(keys) != len(set(keys)):
        raise ValueError("Direct C++ map protocol contains duplicate full identities")
    primary = [item for item in resolved if item["parameter_role"] == "primary"]
    sensitivity = [item for item in resolved if item["parameter_role"] == "sensitivity"]
    if len(primary) != expected_primary_count or len(sensitivity) != expected_sensitivity_count:
        raise ValueError(
            "Direct C++ map protocol has unexpected counts: "
            f"primary={len(primary)} sensitivity={len(sensitivity)}"
        )
    orders = {1, 2, 3, 4, 5, 6}
    section45 = [item for item in primary if item["section"] == "4.5"]
    section46 = [item for item in primary if item["section"] == "4.6"]
    if len(section45) != 54 or len(section46) != 36:
        raise ValueError("Direct protocol must contain 54 Section 4.5 and 36 Section 4.6 primary maps")
    groups: dict[tuple[str, str, str, str], set[int]] = {}
    for item in section45:
        groups.setdefault((item["dataset_id"], _window_label(item["window_size"]), item["weighting_id"], item["section"]), set()).add(item["order"])
    if len(groups) != 9 or any(value != orders for value in groups.values()):
        raise ValueError("Each direct Section 4.5 case must contain uniform orders 1 through 6")
    groups = {}
    for item in section46:
        groups.setdefault((item["dataset_id"], _window_label(item["window_size"]), item["weighting_id"], item["section"]), set()).add(item["order"])
    if len(groups) != 6 or any(value != orders for value in groups.values()):
        raise ValueError("Each direct Section 4.6 case must contain nonuniform orders 1 through 6")
    sensitivity_groups: dict[tuple[str, str, str], set[int]] = {}
    for item in sensitivity:
        sensitivity_groups.setdefault((item["dataset_id"], _window_label(item["window_size"]), item["parameter_set"]), set()).add(item["order"])
    if len(sensitivity_groups) != 6 or any(value != {2, 6} for value in sensitivity_groups.values()):
        raise ValueError("Each direct sensitivity case must contain orders 2 and 6 exactly")
    return {
        "map_count": len(resolved),
        "primary_count": len(primary),
        "sensitivity_count": len(sensitivity),
        "section45_primary_count": len(section45),
        "section46_primary_count": len(section46),
        "dataset_count": len({item["dataset_id"] for item in resolved}),
    }


def validate_cpp_feature_map_protocol(
    specs: Sequence[CppFeatureMapSpec],
    *,
    expected_primary_count: int | None = None,
    expected_sensitivity_count: int | None = None,
) -> dict[str, int]:
    """Validate a frozen map matrix before a formal C++ production run."""

    resolved = [dict(_resolve_feature_spec(spec), order=int(spec.order)) for spec in specs]
    full_keys = [
        (
            item["dataset_id"],
            _window_label(item["window_size"]),
            item["weighting_id"],
            item["parameter_set"],
            item["parameter_role"],
            item["order"],
        )
        for item in resolved
    ]
    if len(set(full_keys)) != len(full_keys):
        raise ValueError("C++ feature-map protocol contains duplicate map identities")
    primary = [item for item in resolved if item["parameter_role"] == "primary"]
    sensitivity = [item for item in resolved if item["parameter_role"] == "sensitivity"]
    primary_modes = {item["mode"] for item in primary}
    if primary and primary_modes != FORMAL_WEIGHT_MODES:
        raise ValueError(f"Primary C++ maps must contain exactly {sorted(FORMAL_WEIGHT_MODES)}, got {sorted(primary_modes)}")
    if expected_primary_count is not None and len(primary) != int(expected_primary_count):
        raise ValueError(f"Expected {expected_primary_count} primary C++ maps, found {len(primary)}")
    if expected_sensitivity_count is not None and len(sensitivity) != int(expected_sensitivity_count):
        raise ValueError(f"Expected {expected_sensitivity_count} sensitivity C++ maps, found {len(sensitivity)}")
    if expected_primary_count == 162 and expected_sensitivity_count == 324:
        datasets = {item["dataset_id"] for item in resolved}
        windows = {_window_label(item["window_size"]) for item in resolved}
        if len(datasets) != 3 or len(windows) != 3:
            raise ValueError("The full formal C++ map protocol requires exactly three datasets and three windows")
        expected_orders = {1, 2, 3, 4, 5, 6}
        primary_by_case: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
        for item in primary:
            primary_by_case.setdefault(
                (item["dataset_id"], _window_label(item["window_size"]), item["mode"]), []
            ).append(item)
        for dataset_id in datasets:
            for window in windows:
                for mode in FORMAL_WEIGHT_MODES:
                    group = primary_by_case.get((dataset_id, window, mode), [])
                    if len(group) != 6 or {item["order"] for item in group} != expected_orders:
                        raise ValueError(
                            f"Primary C++ map protocol is incomplete for {dataset_id}/{window}/{mode}; "
                            "expected one map at each order 1--6"
                        )
        sensitivity_sets = {item["parameter_set"] for item in sensitivity}
        if len(sensitivity_sets) != 6:
            raise ValueError(f"The full formal C++ map protocol requires six sensitivity parameter sets, found {len(sensitivity_sets)}")
        sensitivity_by_case: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
        for item in sensitivity:
            sensitivity_by_case.setdefault(
                (item["dataset_id"], _window_label(item["window_size"]), item["parameter_set"]), []
            ).append(item)
        for dataset_id in datasets:
            for window in windows:
                for parameter_set in sensitivity_sets:
                    group = sensitivity_by_case.get((dataset_id, window, parameter_set), [])
                    if len(group) != 6 or {item["order"] for item in group} != expected_orders:
                        raise ValueError(
                            f"Sensitivity C++ map protocol is incomplete for {dataset_id}/{window}/{parameter_set}; "
                            "expected one map at each order 1--6"
                        )
    return {
        "map_count": len(resolved),
        "primary_count": len(primary),
        "sensitivity_count": len(sensitivity),
        "dataset_count": len({item["dataset_id"] for item in resolved}),
        "window_count": len({_window_label(item["window_size"]) for item in resolved}),
    }


def run_cpp_feature_matrix(
    executable: str | Path,
    specs: Sequence[CppFeatureMapSpec],
    run_dir: str | Path,
    *,
    timeout_seconds: float = 600.0,
    manifest_name: str = "cpp_feature_map_manifest.csv",
    require_valid_support: bool = True,
) -> CppFeatureMatrixResult:
    """Generate coordinate-validated C++ cumulant maps for formal Sections 4.5--4.6.

    The result is the sole C++ production-map interface.  It writes the full
    manifest together with deterministic uniform and nonuniform views, and
    returns every map as a coordinate-indexed dataframe paired with the C++
    support sidecar.  Python analyses should treat ``map_frame`` as the
    authoritative cumulant source and use ``support_frame`` for formal support
    membership rather than recomputing either by output row order.
    """

    executable_path = Path(executable).resolve()
    if not executable_path.is_file():
        raise FileNotFoundError(f"movewindow_cli executable does not exist: {executable_path}")
    target_dir = Path(run_dir).resolve()
    target_dir.mkdir(parents=True, exist_ok=True)
    manifest_file = _safe_token(manifest_name, field="manifest_name")
    if not manifest_file.lower().endswith(".csv"):
        raise ValueError("manifest_name must end in .csv")
    if not specs:
        raise ValueError("At least one formal C++ feature-map specification is required")

    resolved_specs: list[tuple[CppFeatureMapSpec, dict[str, Any]]] = []
    seen_keys: set[tuple[str, str, str, str, str, int]] = set()
    for spec in specs:
        resolved = _resolve_feature_spec(spec)
        identity = (
            resolved["dataset_id"],
            _window_label(resolved["window_size"]),
            resolved["weighting_id"],
            resolved["parameter_set"],
            resolved["parameter_role"],
            int(spec.order),
        )
        if identity in seen_keys:
            raise ValueError(f"Duplicate C++ feature-map identity: {identity}")
        seen_keys.add(identity)
        resolved_specs.append((spec, resolved))

    from window_support import extract_complete_windows

    cli_hash = sha256(executable_path)
    artifacts: list[CppFeatureMapArtifact] = []
    manifest_rows: list[dict[str, Any]] = []
    for spec, resolved in resolved_specs:
        case_id = _feature_case_id(resolved, int(spec.order))
        result = run_cli_case(
            executable_path,
            spec.grid,
            resolved["dataset_id"],
            resolved["window_size"],
            resolved["move_step"],
            resolved["route"],
            resolved["weight_config"],
            resolved["source_spacing"],
            target_dir,
            timeout_seconds,
            int(spec.order),
            case_id=case_id,
            target_property=f"kappa{int(spec.order)}",
            require_window_support_sidecar=True,
        )
        if result.sidecar_frame is None or result.sidecar_path is None:
            raise RuntimeError(f"C++ map {case_id} did not expose a support sidecar")
        batch = extract_complete_windows(spec.grid, resolved["window_size"], resolved["move_step"])
        alignment = validate_cli_window_alignment(
            result.frame,
            result.sidecar_frame,
            batch.coordinates,
            require_valid_support=require_valid_support,
        )
        metadata = _read_metadata(result.metadata_path)
        formula_version = metadata.get("weighting_formula_version")
        if not formula_version:
            raise RuntimeError(f"C++ map {case_id} metadata lacks weighting_formula_version")
        centering_version = metadata.get("cumulant_centering_version")
        if centering_version != CUMULANT_CENTERING_VERSION:
            raise RuntimeError(
                f"C++ map {case_id} has unsupported cumulant_centering_version "
                f"{centering_version!r}; expected {CUMULANT_CENTERING_VERSION!r}"
            )
        expected_centering = (
            "none"
            if int(spec.order) == 1
            else ("window_mean" if resolved["route"] != "weighted" else "weighted_window_mean")
        )
        centering = metadata.get("cumulant_centering")
        if centering != expected_centering:
            raise RuntimeError(
                f"C++ map {case_id} reports cumulant_centering={centering!r}; "
                f"expected {expected_centering!r}"
            )
        valid_support = result.sidecar_frame.loc[
            result.sidecar_frame["map_included"] & result.sidecar_frame["valid"]
        ]
        row = {
            "section": resolved["section"],
            "dataset_id": resolved["dataset_id"],
            "window": _window_label(resolved["window_size"]),
            # ``window`` is kept for older diagnostic readers.  The formal
            # evidence schema uses the explicit stable identifier below.
            "window_id": _window_label(resolved["window_size"]),
            "window_size": "x".join(str(value) for value in resolved["window_size"]),
            "move_step": "x".join(str(value) for value in resolved["move_step"]),
            "weighting": resolved["mode"],
            "weighting_id": resolved["weighting_id"],
            "parameter_set": resolved["parameter_set"],
            "parameter_role": resolved["parameter_role"],
            "order": int(spec.order),
            "algorithm": resolved["algorithm"],
            "authoritative_cplusplus": True,
            "map_node_count": int(len(result.frame)),
            "valid_support_node_count": int(len(valid_support)),
            "map_path": str(result.output_path.resolve()),
            "map_sha256": sha256(result.output_path),
            "sidecar_path": str(result.sidecar_path.resolve()),
            "sidecar_sha256": sha256(result.sidecar_path),
            "metadata_path": str(result.metadata_path.resolve()),
            "metadata_sha256": sha256(result.metadata_path),
            "config_path": str(result.config_path.resolve()),
            "config_sha256": sha256(result.config_path),
            "cli_path": str(executable_path),
            "cli_sha256": cli_hash,
            "weighting_formula_version": formula_version,
            "cumulant_centering": centering,
            "cumulant_centering_version": centering_version,
            "coordinate_validation": str(alignment["status"]),
            "coordinate_validation_map_nodes": int(alignment["map_node_count"]),
            "coordinate_validation_expected_nodes": int(alignment["expected_node_count"]),
            "case_id": result.case_id,
            "cli_elapsed_seconds": float(result.elapsed_seconds),
        }
        artifact = CppFeatureMapArtifact(
            spec=spec,
            cli_result=result,
            map_frame=result.frame.copy(),
            support_frame=result.sidecar_frame.copy(),
            algorithm=resolved["algorithm"],
            coordinate_validation=str(alignment["status"]),
            manifest_row=row,
            authoritative_cplusplus=True,
        )
        artifacts.append(artifact)
        manifest_rows.append(row)

    manifest = pd.DataFrame(manifest_rows).sort_values(
        ["section", "dataset_id", "window", "weighting", "parameter_set", "order"]
    ).reset_index(drop=True)
    manifest_path = target_dir / manifest_file
    uniform_manifest_path = target_dir / "uniform_feature_map_manifest.csv"
    weighted_manifest_path = target_dir / "weighted_feature_map_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    manifest.loc[manifest["weighting"] == "uniform"].to_csv(uniform_manifest_path, index=False)
    manifest.loc[manifest["weighting"] != "uniform"].to_csv(weighted_manifest_path, index=False)
    return CppFeatureMatrixResult(
        artifacts=tuple(artifacts),
        manifest=manifest,
        manifest_path=manifest_path,
        uniform_manifest_path=uniform_manifest_path,
        weighted_manifest_path=weighted_manifest_path,
    )


def compare_maps(reference: pd.DataFrame, candidate: pd.DataFrame, atol: float = ATOL, rtol: float = RTOL) -> dict:
    ref, cand = align_maps(reference, candidate)
    absolute = np.abs(cand - ref)
    allowed = atol + rtol * np.abs(ref)
    failed = absolute > allowed
    return {
        "count": int(len(ref)),
        "failed_count": int(np.sum(failed)),
        "max_abs_error": float(np.max(absolute)) if len(absolute) else 0.0,
        "max_relative_error": float(np.max(absolute / np.maximum(np.abs(ref), atol))) if len(absolute) else 0.0,
        "passed": bool(not np.any(failed)),
    }
