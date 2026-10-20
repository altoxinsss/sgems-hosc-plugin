from __future__ import annotations

"""Compatibility helpers for the Section 4.4 directional-template experiment.

The executable entrypoint delegates to :mod:`formal_runner`, which is the
strict signed-integer, exhaustive-template component runner.  The helpers in
this module remain available to tests and the formal runner; its historical
``_run`` implementation is module-only and must not be used to create a
Section 4.4 evidence package.
"""

import argparse
import csv
import json
import math
import platform
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable

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


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class CliResult:
    condition: str
    case_id: str
    output_grid: Path
    output_meta: Path
    support_csv: Path
    stdout: Path
    stderr: Path
    summary: dict[str, Any]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("run_%Y%m%d_%H%M%S")


def _sha256_json(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _resolve(base: Path, text: str) -> Path:
    path = Path(text)
    return path if path.is_absolute() else (base / path).resolve()


def _load_config(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("experiment configuration must be a JSON object")
    return payload


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in materialized:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(materialized)


def _read_volume(path: Path, property_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    grid = read_columnar_grid(path, property_name)
    return grid.values, grid.x, grid.y, grid.z


def _format_ini_value(value: Iterable[int | float]) -> str:
    return ",".join(str(item) for item in value)


def _write_cli_ini(
    path: Path,
    source_grid: Path,
    source_property: str,
    target_property: str,
    output_dir: Path,
    nod: list[int],
    nlag: list[int],
    lag: list[int],
    azimuth: list[int],
    dip: list[int],
    displayed_dirs: list[int],
    coordinate_dimensions: int,
    source_spacing: list[float],
    direction_vectors: list[list[int]] | list[tuple[int, ...]] | None = None,
    algorithm: str = "AVL Based algorithm",
    template_sampling: str = "exhaustive_integer_offsets",
    fixed_positions: list[int] | None = None,
) -> None:
    if not (len(nod) == len(nlag) == len(lag) == len(azimuth) == len(dip)):
        raise ValueError("all directional CLI parameter lists must have the same length")
    count = len(nod)
    if not displayed_dirs or len(displayed_dirs) > 3:
        raise ValueError("template_cli output must display one to three directions")
    if template_sampling != "exhaustive_integer_offsets":
        raise ValueError("the Section 4.4 runner refuses legacy template search")
    if direction_vectors is None or len(direction_vectors) != count:
        raise ValueError("one signed integer direction vector is required per direction")
    if any(len(vector) != coordinate_dimensions for vector in direction_vectors):
        raise ValueError("integer direction vector dimensions must match the source grid")
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "[general]",
        "calculate_target=Cumulant",
        f"calculate_algorithm={algorithm}",
        f"source_grid={source_grid.as_posix()}",
        f"source_property={source_property}",
        f"coordinate_dimensions={coordinate_dimensions}",
        f"target_property={target_property}",
        "target_origin=0,0,0",
        "source_origin=0,0,0",
        "",
        "[template_shape]",
        f"template_sampling={template_sampling}",
        f"NOD={_format_ini_value(nod)}",
        f"NLAG={_format_ini_value(nlag)}",
        f"LAG={_format_ini_value(lag)}",
        f"LTOL={_format_ini_value([0] * count)}",
        f"AZM={_format_ini_value(azimuth)}",
        f"ATOL={_format_ini_value([0] * count)}",
        f"BANDWH={_format_ini_value([0] * count)}",
        f"DIP={_format_ini_value(dip)}",
        f"DTOL={_format_ini_value([0] * count)}",
        f"BANDWD={_format_ini_value([0] * count)}",
        "integer_direction_vectors=" + ";".join(_format_ini_value(vector) for vector in direction_vectors),
        f"source_spacing={_format_ini_value(source_spacing)}",
        "",
        "[target_cross_sections]",
        f"displayed_dirs={_format_ini_value(displayed_dirs)}",
        f"fixed_positions={_format_ini_value(fixed_positions or [])}",
        "",
        "[runtime]",
        f"output_dir={output_dir.as_posix()}",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def _summary_from_cli_grid(path: Path, property_name: str) -> dict[str, float | int]:
    grid = read_columnar_grid(path, property_name)
    values = grid.values
    return {
        "response_mean": float(np.mean(values)),
        "response_median": float(np.median(values)),
        "response_mean_absolute": float(np.mean(np.abs(values))),
        "response_maximum_absolute": float(np.max(np.abs(values))),
        "response_count": int(values.size),
        "output_nx": int(values.shape[2]),
        "output_ny": int(values.shape[1]),
        "output_nz": int(values.shape[0]),
    }


def _run_cli(
    executable: Path,
    ini_path: Path,
    output_dir: Path,
    property_name: str,
    condition: str,
    case_id: str,
    logs_dir: Path,
    timeout_seconds: int,
    run_label: str = "",
) -> CliResult:
    logs_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"__{run_label}" if run_label else ""
    stdout_path = logs_dir / f"{condition}__{case_id}{suffix}.stdout.txt"
    stderr_path = logs_dir / f"{condition}__{case_id}{suffix}.stderr.txt"
    completed = subprocess.run(
        [str(executable), str(ini_path)],
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=False,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"template_cli failed for {condition}/{case_id} with exit code {completed.returncode}; "
            f"see {stderr_path}"
        )
    output_grid = output_dir / f"{property_name}.txt"
    output_meta = output_dir / f"{property_name}.meta.txt"
    support_csv = output_dir / f"{property_name}.support.csv"
    if not output_grid.exists() or not output_meta.exists() or not support_csv.exists():
        raise RuntimeError(f"template_cli reported success but did not write complete output for {condition}/{case_id}")
    return CliResult(
        condition=condition,
        case_id=case_id,
        output_grid=output_grid,
        output_meta=output_meta,
        support_csv=support_csv,
        stdout=stdout_path,
        stderr=stderr_path,
        summary=_summary_from_cli_grid(output_grid, property_name),
    )


def _run_native_checks(config: dict[str, Any], config_dir: Path, logs_dir: Path) -> list[dict[str, Any]]:
    """Run the configured native-core checks and preserve their complete logs."""

    working_directory = _resolve(config_dir, config["native_test_workdir"])
    logs_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for specification in config["native_test_commands"]:
        name = str(specification["name"])
        command = [str(part) for part in specification["command"]]
        command_path = Path(command[0])
        if not command_path.is_absolute():
            command[0] = str((working_directory / command_path).resolve())
        stdout_path = logs_dir / f"native_{name}.stdout.txt"
        stderr_path = logs_dir / f"native_{name}.stderr.txt"
        completed = subprocess.run(command, cwd=working_directory, text=True, capture_output=True, check=False)
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"native core check '{name}' failed with exit code {completed.returncode}; see {stderr_path}")
        results.append({
            "name": name,
            "command": command,
            "working_directory": str(working_directory),
            "returncode": completed.returncode,
            "stdout": stdout_path,
            "stderr": stderr_path,
        })
    return results


def _relative_l2(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left.ravel())), np.finfo(float).eps)
    return float(np.linalg.norm((left - right).ravel()) / scale)


def _integer_direction(azimuth_degrees: int) -> tuple[int, int]:
    """Map a named azimuth to its signed formal integer direction."""

    return signed_integer_direction(int(azimuth_degrees))


def _matching_explicit_surface(field: np.ndarray, case: dict[str, Any], cli_grid: np.ndarray) -> np.ndarray:
    """Reference the configured two-direction NOD=(2,1) CLI surface explicitly.

    This gate deliberately covers a configuration whose directional template is
    expressible as four explicit offsets: origin, direction 1 twice, and
    direction 2 once. It is separate from the line/rectangle/L reference
    templates reported elsewhere in the experiment.
    """

    if list(case["nod"]) != [2, 1] or len(case["azimuth"]) != 2 or len(case["lag_cells"]) != 2:
        raise ValueError("CLI/reference agreement requires a two-direction NOD=(2,1) case")
    direction_one = _integer_direction(int(case["azimuth"][0]))
    direction_two = _integer_direction(int(case["azimuth"][1]))
    lag_one, lag_two = (int(value) for value in case["lag_cells"])
    response = np.empty(cli_grid.shape, dtype=float)
    for iy in range(cli_grid.shape[0]):
        for ix in range(cli_grid.shape[1]):
            offset_one = (ix * lag_one * direction_one[0], ix * lag_one * direction_one[1])
            offset_two = (iy * lag_two * direction_two[0], iy * lag_two * direction_two[1])
            offsets = np.asarray(((0, 0), offset_one, offset_one, offset_two), dtype=int)
            response[iy, ix] = template_response(field, offsets)
    return response


def _nominal_physical_lag(azimuth_degrees: int, dip_degrees: int, lag_cells: int, spacing: list[float]) -> dict[str, float]:
    """Record a nominal physical lag without claiming a vendor search convention."""

    azimuth = np.deg2rad(float(azimuth_degrees))
    dip = np.deg2rad(float(dip_degrees))
    direction = np.array((np.cos(dip) * np.cos(azimuth), np.cos(dip) * np.sin(azimuth), np.sin(dip)))
    physical = direction * np.asarray(spacing, dtype=float) * float(lag_cells)
    return {
        "nominal_dx": float(physical[0]),
        "nominal_dy": float(physical[1]),
        "nominal_dz": float(physical[2]),
        "nominal_distance": float(np.linalg.norm(physical)),
    }


def _write_physical_lag_metadata(path: Path, config: dict[str, Any]) -> list[dict[str, Any]]:
    spacing = [float(value) for value in config["source_spacing"]]
    rows: list[dict[str, Any]] = []
    for section, cases in (("two_dimensional", config["directional_cases"]), ("three_dimensional", config["three_dimensional_cases"])):
        for case in cases:
            for index, (azimuth, dip, lag) in enumerate(zip(case["azimuth"], case["dip"], case["lag_cells"]), start=1):
                rows.append({
                    "section": section,
                    "case_id": case["id"],
                    "direction_index": index,
                    "azimuth_degrees": azimuth,
                    "dip_degrees": dip,
                    "lag_cells": lag,
                    **_nominal_physical_lag(int(azimuth), int(dip), int(lag), spacing),
                })
    _write_csv(path, rows)
    return rows


def _template_offsets(config: dict[str, Any]) -> dict[str, np.ndarray]:
    output: dict[str, np.ndarray] = {}
    for name, offsets in config["explicit_templates_2d"].items():
        array = np.asarray(offsets, dtype=int)
        if array.shape != (4, 2):
            raise ValueError(f"explicit 2D template '{name}' must contain four x/y offsets")
        output[name] = array
    return output


def _template_offsets_3d(config: dict[str, Any]) -> dict[str, np.ndarray]:
    output: dict[str, np.ndarray] = {}
    for name, offsets in config["explicit_templates_3d"].items():
        array = np.asarray(offsets, dtype=int)
        if array.shape != (4, 3):
            raise ValueError(f"explicit 3D template '{name}' must contain four x/y/z offsets")
        output[name] = array
    return output


def _explicit_2d_rows(
    dataset_id: str,
    condition: str,
    field: np.ndarray,
    templates: dict[str, np.ndarray],
    directions: list[int],
    lags: list[int],
    source_kind: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for template_name, offsets in templates.items():
        for angle in directions:
            for lag in lags:
                rotated = transform_offsets_2d(offsets, angle, lag)
                rows.append({
                    "dataset_id": dataset_id,
                    "condition": condition,
                    "source_kind": source_kind,
                    "template_id": template_name,
                    "angle_degrees": angle,
                    "lag_cells": lag,
                    "calculation_path": "explicit_offset_reference",
                    "response": template_response(field, rotated),
                })
    return rows


def _explicit_3d_rows(dataset_id: str, condition: str, volume: np.ndarray, templates: dict[str, np.ndarray], source_kind: str) -> list[dict[str, Any]]:
    return [
        {
            "dataset_id": dataset_id,
            "condition": condition,
            "source_kind": source_kind,
            "template_id": name,
            "calculation_path": "explicit_offset_reference",
            "response": template_response_nd(volume, offsets),
        }
        for name, offsets in templates.items()
    ]


def _render_figures(run_dir: Path, cli_rows: list[dict[str, Any]], explicit_rows: list[dict[str, Any]],
                    connectivity_rows: list[dict[str, Any]], three_d_cli: list[tuple[str, CliResult]]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures_dir = run_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    output_paths: list[Path] = []

    dataset_ids = sorted({str(row["dataset_id"]) for row in cli_rows})
    for dataset_id in dataset_ids:
        original_rows = [
            row for row in cli_rows
            if row["dataset_id"] == dataset_id and row["source_kind"] == "original"
        ]
        if original_rows:
            fig, axis = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
            for row in original_rows:
                axis.scatter(row["orientation_label"], row["response_mean_absolute"], s=32, label=row["case_id"])
            axis.set_ylabel("Mean absolute fourth-order response")
            axis.set_xlabel("Directional lag-domain configuration")
            axis.tick_params(axis="x", rotation=25)
            axis.legend(fontsize=8)
            path = figures_dir / f"{dataset_id}__directional_cli_response_summary.png"
            fig.savefig(path, dpi=320)
            plt.close(fig)
            output_paths.append(path)

        originals = [
            row for row in explicit_rows
            if row["dataset_id"] == dataset_id and row["source_kind"] == "original"
        ]
        surrogates = [
            row for row in explicit_rows
            if row["dataset_id"] == dataset_id and row["source_kind"] == "iaaft_surrogate"
        ]
        if originals and surrogates:
            figure, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)
            for axis, template_name in zip(axes, sorted({row["template_id"] for row in originals})[:2]):
                baseline = [row for row in originals if row["template_id"] == template_name]
                null = [row["response"] for row in surrogates if row["template_id"] == template_name]
                axis.hist(null, bins=min(16, max(6, len(null) // 2)), color="#9ba9b4", edgecolor="black")
                if baseline:
                    axis.axvline(baseline[0]["response"], color="#b62f28", linewidth=1.6)
                axis.set_title(template_name)
                axis.set_xlabel("Fourth-order response")
                axis.set_ylabel("IAAFT fields")
            path = figures_dir / f"{dataset_id}__explicit_template_iaaft_distribution.png"
            figure.savefig(path, dpi=320)
            plt.close(figure)
            output_paths.append(path)

        original_connectivity = [
            row for row in connectivity_rows
            if row["dataset_id"] == dataset_id and row["source_kind"] == "original" and row["dimension"] == "2d"
        ]
        surrogate_connectivity = [
            row for row in connectivity_rows
            if row["dataset_id"] == dataset_id and row["source_kind"] == "iaaft_surrogate" and row["dimension"] == "2d"
        ]
        if original_connectivity and surrogate_connectivity:
            figure, axis = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
            labels = ["original", "IAAFT mean"]
            values = [
                original_connectivity[0]["largest_component_fraction_of_active"],
                float(np.mean([row["largest_component_fraction_of_active"] for row in surrogate_connectivity])),
            ]
            axis.bar(labels, values, color=("#3274a1", "#a6a6a6"), edgecolor="black")
            axis.set_ylim(0.0, 1.0)
            axis.set_ylabel("Largest component / active cells")
            path = figures_dir / f"{dataset_id}__connectivity_largest_component.png"
            figure.savefig(path, dpi=320)
            plt.close(figure)
            output_paths.append(path)

    for dataset_id, result in three_d_cli:
        grid = read_columnar_grid(result.output_grid)
        middle = grid.values.shape[0] // 2
        figure, axis = plt.subplots(figsize=(5.0, 4.2), constrained_layout=True)
        image = axis.imshow(grid.values[middle], origin="lower", cmap="coolwarm")
        axis.set_xlabel("Lag-direction 1 index")
        axis.set_ylabel("Lag-direction 2 index")
        axis.set_title(f"{dataset_id}: 3D response at central lag-direction 3")
        figure.colorbar(image, ax=axis, label="Fourth-order response")
        path = figures_dir / f"{dataset_id}__three_dimensional_template_response_slice.png"
        figure.savefig(path, dpi=320)
        plt.close(figure)
        output_paths.append(path)
    return output_paths


def _module_only_status(gates: dict[str, bool]) -> str:
    """Return a non-final status for the retired module-only execution path."""

    return "DIAGNOSTIC" if all(gates.values()) else "FAIL"


def _write_module_report(run_dir: Path, status: str, gates: dict[str, bool], config: dict[str, Any], cli_rows: list[dict[str, Any]],
                         surrogate_audit: list[dict[str, Any]], rotation_rows: list[dict[str, Any]]) -> Path:
    report = run_dir / "MODULE_RUN_REPORT.md"
    accepted = sum(1 for row in surrogate_audit if bool(row["accepted"]))
    lines = [
        "# Section 4.4 Directional Template Module Diagnostic",
        "",
        f"Execution status: **{status}**",
        "",
        "## Gates",
        "",
    ]
    lines.extend(f"- `{name}`: `{'PASS' if value else 'FAIL'}`" for name, value in gates.items())
    lines.extend([
        "",
        "## Output Layout",
        "",
        "- `configs/`: generated production `template_cli` configurations.",
        "- `inputs/`: original, rotated, and accepted IAAFT grid inputs used by CLI.",
        "- `cli_outputs/`: direct production CLI output grids and metadata.",
        "- `tables/`: numerical response, rotation, surrogate-gate, connectivity, and nominal-lag tables.",
        "- `figures/`: publication-oriented diagnostic figures generated from the tables.",
        "- `evidence_manifest.json`: source, executable, configuration, environment, and output hashes.",
        "",
        "## Interpretation Boundary",
        "",
        "Rows labelled `production_template_cli` were calculated by the production CLI. Rows labelled "
        "`explicit_offset_reference` use recorded offsets because the CLI accepts directional lag configurations rather "
        "than arbitrary explicit point lists. IAAFT acceptance is determined solely by histogram and variogram gates.",
        "",
        "## Counts",
        "",
        f"- Production CLI response summaries: {len(cli_rows)}",
        f"- IAAFT attempts: {len(surrogate_audit)}",
        f"- Accepted IAAFT fields: {accepted}",
        f"- Rotation comparisons: {len(rotation_rows)}",
        f"- Requested IAAFT fields per input: {config['iaaft']['count']}",
    ])
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def _write_module_results_index(run_dir: Path) -> Path:
    path = run_dir / "MODULE_RESULTS_INDEX.md"
    path.write_text(
        "# Section 4.4 Module Diagnostic Index\n\n"
        "This module-only output is not a formal Section 4.4 evidence package. "
        "Only the top-level formal entrypoint can create the study result index and final PASS status.\n\n"
        "## 4.4.1 Two-dimensional directional and distance response\n\n"
        "Use `tables/production_cli_response_summary.csv` and "
        "`figures/*__directional_cli_response_summary.png` for lag-domain response magnitude. "
        "Use `tables/explicit_template_responses.csv` for line, rectangle, and L-template offset responses.\n\n"
        "## 4.4.2 Rotation, surrogate, and connectivity checks\n\n"
        "Use `tables/rotation_validation.csv` for original versus 90-degree rotated directional-lag surfaces. "
        "Use `tables/iaaft_acceptance.csv`, `figures/*__explicit_template_iaaft_distribution.png`, and "
        "`tables/connectivity_metrics.csv` to distinguish histogram/variogram preservation from connectivity and "
        "higher-order template response. Use `tables/cli_reference_agreement.csv` for the matching production CLI "
        "and explicit-offset directional surface gate.\n\n"
        "## 4.4.3 Three-dimensional templates and physical lag metadata\n\n"
        "Use `tables/nominal_physical_lags.csv`, `tables/three_dimensional_cli_response_summary.csv`, and "
        "`figures/*__three_dimensional_template_response_slice.png`. The nominal physical lag table records configuration "
        "angles, cell lags, spacing, and the derived nominal vector; it does not replace the CLI search metadata.\n\n"
        "## Reproducibility\n\n"
        "Use `MODULE_RUN_REPORT.md`, `evidence_manifest.json`, generated `configs/`, and `logs/` to trace this diagnostic run.\n",
        encoding="utf-8",
    )
    return path


def _collect_hashes(root: Path) -> dict[str, str]:
    return {
        file.relative_to(root).as_posix(): sha256_file(file)
        for file in sorted(root.rglob("*"))
        if file.is_file() and file.name != "evidence_manifest.json"
    }


def _run(config_path: Path, output_root: Path | None, smoke: bool) -> Path:
    config = _load_config(config_path)
    config_dir = config_path.parent
    executable = _resolve(config_dir, config["template_cli"])
    if not executable.exists():
        raise FileNotFoundError(f"template_cli does not exist: {executable}")
    datasets = list(config["datasets"])
    directional_cases = list(config["directional_cases"])
    three_d_cases = list(config["three_dimensional_cases"])
    rotation_pairs_config = list(config["rotation_validation"]["pairs"])
    iaaft = dict(config["iaaft"])
    if smoke:
        datasets = datasets[:1]
        # Keep every endpoint required by the configured compact rotation
        # checks. This avoids a smoke run silently dropping a comparison.
        selected_case_ids = {
            str(pair[key])
            for pair in rotation_pairs_config
            for key in ("original_case_id", "rotated_case_id")
            if bool(pair.get("include_in_smoke", False))
        }
        directional_cases = [case for case in directional_cases if str(case["id"]) in selected_case_ids]
        three_d_cases = three_d_cases[:1]
        iaaft["count"] = min(2, int(iaaft["count"]))
        iaaft["max_attempts"] = min(12, int(iaaft["max_attempts"]))
        iaaft["iterations"] = min(40, int(iaaft["iterations"]))

    root = output_root or _resolve(config_dir, config["output_root"])
    run_dir = root / _timestamp()
    for name in ("configs", "inputs", "cli_outputs", "logs", "tables", "figures"):
        (run_dir / name).mkdir(parents=True, exist_ok=False)
    shutil.copy2(config_path, run_dir / "experiment_config.json")
    native_checks = _run_native_checks(config, config_dir, run_dir / "logs")
    preprocessing_manifest = _resolve(config_dir, config["preprocessing_manifest"])
    if not preprocessing_manifest.exists():
        raise FileNotFoundError(f"preprocessing manifest does not exist: {preprocessing_manifest}")

    source_spacing = [float(value) for value in config["source_spacing"]]
    templates_2d = _template_offsets(config)
    templates_3d = _template_offsets_3d(config)
    _write_physical_lag_metadata(run_dir / "tables" / "nominal_physical_lags.csv", {
        **config, "directional_cases": directional_cases, "three_dimensional_cases": three_d_cases,
    })

    cli_rows: list[dict[str, Any]] = []
    explicit_rows: list[dict[str, Any]] = []
    connectivity_rows: list[dict[str, Any]] = []
    surrogate_audit: list[dict[str, Any]] = []
    three_d_cli: list[tuple[str, CliResult]] = []
    input_hashes: dict[str, str] = {}
    two_d_fields: dict[str, np.ndarray] = {}

    for dataset_index, dataset in enumerate(datasets):
        grid_path = _resolve(config_dir, dataset["grid"])
        volume, x, y, z = _read_volume(grid_path, dataset["property"])
        if volume.shape[0] < 4:
            raise ValueError(f"3D dataset '{dataset['id']}' requires at least four z layers")
        target_z = int(dataset["slice_z"])
        matching = np.flatnonzero(z == target_z)
        if matching.size != 1:
            raise ValueError(f"slice z={target_z} is not present in dataset '{dataset['id']}'")
        field = volume[int(matching[0])]
        two_d_fields[str(dataset["id"])] = field
        # Array rows increase with the grid y coordinate. A clockwise array
        # rotation is therefore a +90-degree rotation in the x/y convention
        # used by template_cli: +x maps to +y and +y maps to -x.
        rotated = np.rot90(field, k=-1)
        source_prefix = dataset["id"]
        original_path = run_dir / "inputs" / f"{source_prefix}__original.grid"
        rotated_path = run_dir / "inputs" / f"{source_prefix}__rotated90.grid"
        write_columnar_grid(original_path, field, dataset["property"])
        write_columnar_grid(rotated_path, rotated, dataset["property"])
        input_hashes[grid_path.as_posix()] = sha256_file(grid_path)
        input_hashes[original_path.relative_to(run_dir).as_posix()] = sha256_file(original_path)
        input_hashes[rotated_path.relative_to(run_dir).as_posix()] = sha256_file(rotated_path)

        originals = [("original", original_path, field), ("rotated90", rotated_path, rotated)]
        for condition, source_path, condition_field in originals:
            explicit_rows.extend(_explicit_2d_rows(
                source_prefix, f"{source_prefix}__{condition}", condition_field, templates_2d,
                config["explicit_directions_degrees"], config["explicit_lags_cells"], condition,
            ))
            mask = condition_field >= float(config["connectivity_threshold"])
            connectivity_rows.append({
                "dataset_id": source_prefix, "condition": condition, "source_kind": condition,
                "dimension": "2d", **connectivity_metrics(mask, ("y", "x")),
            })
            results_by_case: dict[str, CliResult] = {}
            for case in directional_cases:
                case_id = str(case["id"])
                property_name = f"{source_prefix}__{condition}__{case_id}"
                output_dir = run_dir / "cli_outputs" / source_prefix / condition / case_id
                ini_path = run_dir / "configs" / f"{property_name}.ini"
                _write_cli_ini(
                    ini_path, source_path, dataset["property"], property_name, output_dir,
                    [int(value) for value in case["nod"]], [int(value) for value in case["nlag"]],
                    [int(value) for value in case["lag_cells"]], [int(value) for value in case["azimuth"]],
                    [int(value) for value in case["dip"]], [int(value) for value in case["displayed_dirs"]],
                    2, source_spacing,
                )
                result = _run_cli(executable, ini_path, output_dir, property_name, condition, case_id,
                                  run_dir / "logs" / source_prefix, int(config["template_cli_timeout_seconds"]))
                results_by_case[case_id] = result
                cli_rows.append({
                    "dataset_id": source_prefix, "condition": condition, "source_kind": condition,
                    "case_id": case_id, "orientation_label": case["orientation_label"],
                    "calculation_path": "production_template_cli", "output_grid": result.output_grid.relative_to(run_dir).as_posix(),
                    **result.summary,
                })
            if condition == "original":
                for case in directional_cases:
                    rotated_case_id = str(case["rotated_case_id"])
                    # The paired rotated result is created in the next loop condition; retain source results by ID.
                    results_by_case[case["id"]] = results_by_case[str(case["id"])]

        # Run each accepted IAAFT field only after all gates are evaluated. No response values enter this selection.
        surrogates, audit = generate_accepted_iaaft_surrogates(
            field, int(iaaft["count"]), int(iaaft["seed"]) + dataset_index * 100000,
            int(iaaft["iterations"]), [int(value) for value in iaaft["variogram_lags_cells"]],
            float(iaaft["max_variogram_relative_rmse"]), int(iaaft["max_attempts"]),
        )
        for row in audit:
            surrogate_audit.append({"dataset_id": source_prefix, **row})
        for surrogate_index, surrogate in enumerate(surrogates, start=1):
            condition = f"iaaft_{surrogate_index:02d}"
            source_path = run_dir / "inputs" / f"{source_prefix}__{condition}.grid"
            write_columnar_grid(source_path, surrogate, dataset["property"])
            input_hashes[source_path.relative_to(run_dir).as_posix()] = sha256_file(source_path)
            explicit_rows.extend(_explicit_2d_rows(
                source_prefix, f"{source_prefix}__{condition}", surrogate, templates_2d,
                config["explicit_directions_degrees"], config["explicit_lags_cells"], "iaaft_surrogate",
            ))
            connectivity_rows.append({
                "dataset_id": source_prefix, "condition": condition, "source_kind": "iaaft_surrogate",
                "dimension": "2d", **connectivity_metrics(surrogate >= float(config["connectivity_threshold"]), ("y", "x")),
            })
            for case in directional_cases:
                case_id = str(case["id"])
                property_name = f"{source_prefix}__{condition}__{case_id}"
                output_dir = run_dir / "cli_outputs" / source_prefix / condition / case_id
                ini_path = run_dir / "configs" / f"{property_name}.ini"
                _write_cli_ini(
                    ini_path, source_path, dataset["property"], property_name, output_dir,
                    [int(value) for value in case["nod"]], [int(value) for value in case["nlag"]],
                    [int(value) for value in case["lag_cells"]], [int(value) for value in case["azimuth"]],
                    [int(value) for value in case["dip"]], [int(value) for value in case["displayed_dirs"]],
                    2, source_spacing,
                )
                result = _run_cli(executable, ini_path, output_dir, property_name, "iaaft_surrogate", case_id,
                                  run_dir / "logs" / source_prefix, int(config["template_cli_timeout_seconds"]))
                cli_rows.append({
                    "dataset_id": source_prefix, "condition": condition, "source_kind": "iaaft_surrogate",
                    "case_id": case_id, "orientation_label": case["orientation_label"],
                    "calculation_path": "production_template_cli", "output_grid": result.output_grid.relative_to(run_dir).as_posix(),
                    **result.summary,
                })

        # 3D CLI and explicit templates are evaluated on the unrotated volume.
        explicit_rows.extend(_explicit_3d_rows(source_prefix, source_prefix, volume, templates_3d, "original_3d"))
        connectivity_rows.append({
            "dataset_id": source_prefix, "condition": "original_3d", "source_kind": "original",
            "dimension": "3d", **connectivity_metrics(volume >= float(config["connectivity_threshold"]), ("z", "y", "x")),
        })
        for case in three_d_cases:
            case_id = str(case["id"])
            property_name = f"{source_prefix}__3d__{case_id}"
            output_dir = run_dir / "cli_outputs" / source_prefix / "three_dimensional" / case_id
            ini_path = run_dir / "configs" / f"{property_name}.ini"
            _write_cli_ini(
                ini_path, grid_path, dataset["property"], property_name, output_dir,
                [int(value) for value in case["nod"]], [int(value) for value in case["nlag"]],
                [int(value) for value in case["lag_cells"]], [int(value) for value in case["azimuth"]],
                [int(value) for value in case["dip"]], [int(value) for value in case["displayed_dirs"]],
                3, source_spacing,
            )
            result = _run_cli(executable, ini_path, output_dir, property_name, "original_3d", case_id,
                              run_dir / "logs" / source_prefix, int(config["template_cli_timeout_seconds"]))
            three_d_cli.append((source_prefix, result))
            row = {
                "dataset_id": source_prefix, "condition": "original_3d", "source_kind": "original",
                "case_id": case_id, "orientation_label": case["orientation_label"],
                "calculation_path": "production_template_cli", "output_grid": result.output_grid.relative_to(run_dir).as_posix(),
                **result.summary,
            }
            cli_rows.append(row)

    rotation_rows: list[dict[str, Any]] = []
    for dataset in datasets:
        dataset_id = str(dataset["id"])
        original_by_case = {row["case_id"]: row for row in cli_rows if row["dataset_id"] == dataset_id and row["condition"] == "original"}
        rotated_by_case = {row["case_id"]: row for row in cli_rows if row["dataset_id"] == dataset_id and row["condition"] == "rotated90"}
        for pair in rotation_pairs_config:
            if smoke and not bool(pair.get("include_in_smoke", False)):
                continue
            source_case, target_case = str(pair["original_case_id"]), str(pair["rotated_case_id"])
            if source_case not in original_by_case or target_case not in rotated_by_case:
                raise RuntimeError(f"missing production rotation pair {source_case}->{target_case}")
            left = read_columnar_grid(run_dir / original_by_case[source_case]["output_grid"]).values
            right = read_columnar_grid(run_dir / rotated_by_case[target_case]["output_grid"]).values
            rotation_rows.append({
                "dataset_id": dataset_id,
                "original_case_id": source_case,
                "rotated_case_id": target_case,
                "relative_l2_error": _relative_l2(left, right),
                "maximum_absolute_error": float(np.max(np.abs(left - right))) if left.shape == right.shape else float("inf"),
                "shape_match": bool(left.shape == right.shape),
            })

    agreement_config = dict(config["cli_reference_agreement"])
    agreement_case_id = str(agreement_config["case_id"])
    case_by_id = {str(case["id"]): case for case in directional_cases}
    if agreement_case_id not in case_by_id:
        raise RuntimeError(f"CLI/reference agreement case '{agreement_case_id}' was not selected")
    agreement_rows: list[dict[str, Any]] = []
    agreement_surfaces = run_dir / "tables" / "cli_reference_surfaces"
    agreement_surfaces.mkdir(exist_ok=True)
    for dataset in datasets:
        dataset_id = str(dataset["id"])
        matches = [
            row for row in cli_rows
            if row["dataset_id"] == dataset_id and row["condition"] == "original" and row["case_id"] == agreement_case_id
        ]
        if len(matches) != 1:
            raise RuntimeError(f"expected one original CLI output for {dataset_id}/{agreement_case_id}")
        cli_grid = read_columnar_grid(run_dir / matches[0]["output_grid"]).values[0]
        reference_grid = _matching_explicit_surface(two_d_fields[dataset_id], case_by_id[agreement_case_id], cli_grid)
        surface_path = agreement_surfaces / f"{dataset_id}__{agreement_case_id}.npz"
        np.savez_compressed(surface_path, production_template_cli=cli_grid, explicit_reference=reference_grid)
        agreement_rows.append({
            "dataset_id": dataset_id,
            "case_id": agreement_case_id,
            "production_path": matches[0]["output_grid"],
            "reference_path": surface_path.relative_to(run_dir).as_posix(),
            "response_count": int(cli_grid.size),
            "relative_l2_error": _relative_l2(cli_grid, reference_grid),
            "maximum_absolute_error": float(np.max(np.abs(cli_grid - reference_grid))),
            "max_relative_l2_error": float(agreement_config["max_relative_l2_error"]),
            "max_absolute_error": float(agreement_config["max_absolute_error"]),
        })

    _write_csv(run_dir / "tables" / "production_cli_response_summary.csv", cli_rows)
    _write_csv(run_dir / "tables" / "three_dimensional_cli_response_summary.csv", [
        row for row in cli_rows if row["condition"] == "original_3d"
    ])
    _write_csv(run_dir / "tables" / "explicit_template_responses.csv", explicit_rows)
    _write_csv(run_dir / "tables" / "connectivity_metrics.csv", connectivity_rows)
    _write_csv(run_dir / "tables" / "iaaft_acceptance.csv", surrogate_audit)
    _write_csv(run_dir / "tables" / "rotation_validation.csv", rotation_rows)
    _write_csv(run_dir / "tables" / "cli_reference_agreement.csv", agreement_rows)

    figures = _render_figures(run_dir, cli_rows, explicit_rows, connectivity_rows, three_d_cli)
    max_rotation_error = float(config["rotation_validation"]["max_relative_l2_error"])
    expected_surrogates = len(datasets) * int(iaaft["count"])
    gates = {
        "production_cli_outputs_complete": bool(cli_rows),
        "native_core_checks": bool(native_checks) and all(check["returncode"] == 0 for check in native_checks),
        "iaaft_accepted_count": sum(1 for row in surrogate_audit if bool(row["accepted"])) == expected_surrogates,
        "iaaft_histogram_gate": all(bool(row["histogram_exact"]) for row in surrogate_audit),
        # Rejected IAAFT attempts remain in the audit table. The gate applies
        # to the accepted fields that enter response calculations.
        "iaaft_variogram_gate": all(bool(row["variogram_pass"]) for row in surrogate_audit if bool(row["accepted"])),
        "rotation_surface_shape_gate": all(bool(row["shape_match"]) for row in rotation_rows),
        "rotation_surface_error_gate": all(float(row["relative_l2_error"]) <= max_rotation_error for row in rotation_rows),
        "cli_reference_agreement_gate": all(
            float(row["relative_l2_error"]) <= float(row["max_relative_l2_error"])
            and float(row["maximum_absolute_error"]) <= float(row["max_absolute_error"])
            for row in agreement_rows
        ),
        "three_dimensional_cli_outputs_complete": bool(three_d_cli),
        "figures_complete": len(figures) >= 3,
    }
    status = _module_only_status(gates)
    report = _write_module_report(run_dir, status, gates, {**config, "iaaft": iaaft}, cli_rows, surrogate_audit, rotation_rows)
    results_index = _write_module_results_index(run_dir)
    manifest = {
        "experiment": "template_directional_20261010",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "overall_gate": status,
        "mode": "module_only",
        "gates": gates,
        "config": {"path": str(config_path), "sha256": _sha256_json(config_path)},
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__).resolve())},
        "template_cli": {"path": str(executable), "sha256": sha256_file(executable)},
        "native_core_checks": [
            {
                "name": check["name"],
                "command": check["command"],
                "working_directory": check["working_directory"],
                "stdout": check["stdout"].relative_to(run_dir).as_posix(),
                "stderr": check["stderr"].relative_to(run_dir).as_posix(),
            }
            for check in native_checks
        ],
        "preprocessing_manifest": {"path": str(preprocessing_manifest), "sha256": sha256_file(preprocessing_manifest)},
        "inputs": input_hashes,
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
        "run_parameters": {"smoke": smoke, "iaaft": iaaft, "datasets": [dataset["id"] for dataset in datasets]},
        "reports": [report.relative_to(run_dir).as_posix(), results_index.relative_to(run_dir).as_posix()],
        "outputs": _collect_hashes(run_dir),
    }
    (run_dir / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(run_dir)
    return run_dir


def main() -> int:
    # Keep the historical filename as a stable command path while making the
    # strict integer-offset runner the only executable Section 4.4 entry.
    from formal_runner import main as formal_main

    return formal_main()


if __name__ == "__main__":
    raise SystemExit(main())
