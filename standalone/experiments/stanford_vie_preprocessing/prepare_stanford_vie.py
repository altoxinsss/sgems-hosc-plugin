from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


THIS_DIR = Path(__file__).resolve().parent
R1_5_DIR = THIS_DIR.parent / "r1_5_incremental_value_20260924"
if str(R1_5_DIR) not in sys.path:
    sys.path.insert(0, str(R1_5_DIR))

from grid_io import write_map
from validation import block_reduce_mean, channel_count_2x2x2, read_facies_dat, strict_majority_channel_labels


PROPERTY_NAMES = (
    "porosity",
    "channel_count",
    "channel_occupancy",
    "channel_strict_majority",
    "channel_label_valid",
)


@dataclass(frozen=True)
class CoarseVolumes:
    """Co-located coarse Stanford VI-E arrays in canonical (z, y, x) order."""

    porosity: np.ndarray
    channel_count: np.ndarray
    channel_occupancy: np.ndarray
    channel_strict_majority: np.ndarray
    channel_label_valid: np.ndarray

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(length) for length in self.porosity.shape)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_porosity(path: Path, shape_zyx: tuple[int, int, int], header_lines: int) -> np.ndarray:
    values = np.loadtxt(path, skiprows=header_lines, dtype=float)
    if values.size != int(np.prod(shape_zyx)):
        raise ValueError(f"Porosity value count {values.size} does not match {shape_zyx}")
    return values.reshape(shape_zyx)


def read_poro_inc(path: Path) -> np.ndarray:
    text = path.read_text(encoding="utf-8").replace("PORO", " ").replace("/", " ")
    return np.fromstring(text, sep=" ", dtype=float)


def coarse_volumes(config: dict[str, Any]) -> CoarseVolumes:
    """Apply the frozen transform to porosity and facies inputs.

    Fine input arrays are interpreted as (z, y, x). Both fields are reduced in
    contiguous 2-by-2-by-2 blocks and then reversed along the coarse z and y
    axes. Every emitted property therefore remains co-located by construction.
    """
    fine_shape = tuple(int(value) for value in config["fine_shape_zyx"])
    if len(fine_shape) != 3 or any(length <= 0 or length % 2 for length in fine_shape):
        raise ValueError("fine_shape_zyx must contain three positive even lengths")
    porosity_fine = read_porosity(Path(config["porosity_path"]), fine_shape, int(config["porosity_header_lines"]))
    facies_fine = read_facies_dat(Path(config["facies_path"]), fine_shape, int(config["facies_header_lines"]))
    target_shape = tuple(length // 2 for length in fine_shape)
    porosity = block_reduce_mean(porosity_fine, target_shape)[::-1, ::-1, :]
    channel_count = channel_count_2x2x2(facies_fine, int(config["channel_code"]))[::-1, ::-1, :]
    channel_occupancy = channel_count.astype(float) / 8.0
    labels, valid = strict_majority_channel_labels(channel_count)
    arrays = (porosity, channel_count, channel_occupancy, labels, valid)
    if any(array.shape != target_shape for array in arrays):
        raise RuntimeError("The coarse properties are not co-located")
    if not np.array_equal(channel_occupancy, channel_count / 8.0):
        raise RuntimeError("channel_occupancy is inconsistent with channel_count")
    return CoarseVolumes(porosity, channel_count, channel_occupancy, labels, valid)


def grid_coordinates(z_start: int, z_end: int, shape_yx: tuple[int, int]) -> pd.DataFrame:
    if z_end < z_start:
        raise ValueError("Subvolume end must not precede subvolume start")
    y_length, x_length = (int(value) for value in shape_yx)
    return pd.DataFrame(
        [(x, y, z) for z in range(z_start, z_end + 1) for y in range(1, y_length + 1) for x in range(1, x_length + 1)],
        columns=["x", "y", "z"],
    )


def _subvolume_slice(volumes: CoarseVolumes, z_start: int, z_end: int) -> slice:
    z_length = volumes.shape[0]
    if z_start < 1 or z_end > z_length:
        raise ValueError(f"Subvolume z={z_start}..{z_end} lies outside 1..{z_length}")
    return slice(z_start - 1, z_end)


def _subvolume_arrays(volumes: CoarseVolumes, z_start: int, z_end: int) -> dict[str, np.ndarray]:
    slc = _subvolume_slice(volumes, z_start, z_end)
    return {
        "porosity": volumes.porosity[slc],
        "channel_count": volumes.channel_count[slc],
        "channel_occupancy": volumes.channel_occupancy[slc],
        "channel_strict_majority": volumes.channel_strict_majority[slc].astype(int),
        "channel_label_valid": volumes.channel_label_valid[slc].astype(int),
    }


def write_subvolume(
    output_root: Path,
    volumes: CoarseVolumes,
    z_start: int,
    z_end: int,
) -> tuple[list[Path], dict[str, Any]]:
    """Write one co-located volume and return complete manifest metadata."""
    subdir = output_root / f"z{z_start}_{z_end}"
    subdir.mkdir(parents=True, exist_ok=True)
    arrays = _subvolume_arrays(volumes, z_start, z_end)
    expected_shape = (z_end - z_start + 1, volumes.shape[1], volumes.shape[2])
    if any(array.shape != expected_shape for array in arrays.values()):
        raise RuntimeError("Subvolume property shapes differ")
    coordinates = grid_coordinates(z_start, z_end, (volumes.shape[1], volumes.shape[2]))
    expected_row_count = int(np.prod(expected_shape))
    if len(coordinates) != expected_row_count:
        raise RuntimeError("Coordinate count does not match the subvolume shape")

    paths: list[Path] = []
    output_records: list[dict[str, Any]] = []
    for name in PROPERTY_NAMES:
        values = arrays[name].reshape(-1)
        path = subdir / f"stanford_vie_z{z_start}_{z_end}_{name}.grid"
        write_map(path, coordinates, values, name, title=f"Stanford VI-E z={z_start}..{z_end}; {name}")
        paths.append(path)
        output_records.append(
            {
                "path": str(path.resolve()),
                "property_name": name,
                "row_count": int(values.size),
                "sha256": sha256(path),
            }
        )

    valid = arrays["channel_label_valid"].astype(bool)
    labels = arrays["channel_strict_majority"].astype(bool)
    metadata = {
        "subvolume_id": f"z{z_start}_{z_end}",
        "shape_zyx": list(expected_shape),
        "coordinate_range": {
            "x": [1, volumes.shape[2]],
            "y": [1, volumes.shape[1]],
            "z": [z_start, z_end],
        },
        "row_count": expected_row_count,
        "positive_count": int(np.count_nonzero(labels & valid)),
        "negative_count": int(np.count_nonzero(~labels & valid)),
        "tie_count": int(np.count_nonzero(~valid)),
        "outputs": output_records,
    }
    return paths, metadata


def _fixture_value(volume: np.ndarray, x: int, y: int, z: int) -> Any:
    return volume[z - 1, y - 1, x - 1].item()


def validate_co_location_fixtures(
    volumes: CoarseVolumes,
    fixtures: list[dict[str, Any]],
    *,
    porosity_atol: float,
) -> list[dict[str, Any]]:
    """Validate fixed output coordinates after the Y/Z reversal.

    A fixture includes values from every property at one output coordinate. It
    detects an axis reversal, a property-only transform, or a label alignment
    failure that aggregate counts cannot reveal.
    """
    if not fixtures:
        raise ValueError("At least one co-location fixture is required")
    checked: list[dict[str, Any]] = []
    z_length, y_length, x_length = volumes.shape
    for fixture in fixtures:
        coordinate = fixture.get("coordinate")
        if not isinstance(coordinate, list) or len(coordinate) != 3:
            raise ValueError("Each co-location fixture requires coordinate [x, y, z]")
        x, y, z = (int(value) for value in coordinate)
        if not (1 <= x <= x_length and 1 <= y <= y_length and 1 <= z <= z_length):
            raise ValueError(f"Co-location fixture coordinate {(x, y, z)} is outside the coarse grid")
        observed = {
            "porosity": float(_fixture_value(volumes.porosity, x, y, z)),
            "channel_count": int(_fixture_value(volumes.channel_count, x, y, z)),
            "channel_occupancy": float(_fixture_value(volumes.channel_occupancy, x, y, z)),
            "channel_strict_majority": int(_fixture_value(volumes.channel_strict_majority, x, y, z)),
            "channel_label_valid": int(_fixture_value(volumes.channel_label_valid, x, y, z)),
        }
        expected = fixture.get("expected")
        if not isinstance(expected, dict):
            raise ValueError("Each co-location fixture requires an expected property mapping")
        if set(expected) != set(PROPERTY_NAMES):
            raise ValueError("Each co-location fixture must specify every co-located property")
        if not np.isclose(observed["porosity"], float(expected["porosity"]), rtol=0.0, atol=porosity_atol):
            raise RuntimeError(f"Co-location fixture porosity mismatch at {(x, y, z)}")
        for name in ("channel_count", "channel_strict_majority", "channel_label_valid"):
            if observed[name] != int(expected[name]):
                raise RuntimeError(f"Co-location fixture {name} mismatch at {(x, y, z)}")
        if not np.isclose(
            observed["channel_occupancy"], float(expected["channel_occupancy"]), rtol=0.0, atol=1.0e-15
        ):
            raise RuntimeError(f"Co-location fixture channel_occupancy mismatch at {(x, y, z)}")
        checked.append({"coordinate_xyz": [x, y, z], "status": "PASS", "observed": observed})
    return checked


def load_co_location_fixtures(config: dict[str, Any], config_path: Path) -> tuple[list[dict[str, Any]], Path]:
    raw_path = config.get("co_location_fixture_path")
    if not isinstance(raw_path, str) or not raw_path:
        raise ValueError("co_location_fixture_path is required")
    fixture_path = Path(raw_path)
    if not fixture_path.is_absolute():
        fixture_path = config_path.parent / fixture_path
    payload = json.loads(fixture_path.read_text(encoding="utf-8"))
    fixtures = payload.get("fixtures")
    if not isinstance(fixtures, list):
        raise ValueError("Co-location fixture file requires a fixtures list")
    return fixtures, fixture_path.resolve()


def _axis_mapping(fine_shape: tuple[int, int, int], coarse_shape: tuple[int, int, int]) -> dict[str, Any]:
    return {
        "input_array_order": ["z", "y", "x"],
        "output_array_order": ["z", "y", "x"],
        "fine_shape_zyx": list(fine_shape),
        "coarse_shape_zyx": list(coarse_shape),
        "reduction": "contiguous_2x2x2_blocks; porosity=mean; facies=channel_count",
        "source_reduced_to_output_index": {
            "z": f"output_z_index = {coarse_shape[0] - 1} - source_reduced_z_index",
            "y": f"output_y_index = {coarse_shape[1] - 1} - source_reduced_y_index",
            "x": "output_x_index = source_reduced_x_index",
        },
        "output_coordinate_to_index": {
            "x": "x_index = x_coordinate - 1",
            "y": "y_index = y_coordinate - 1",
            "z": "z_index = z_coordinate - 1",
        },
    }


def build_manifest(
    *,
    config: dict[str, Any],
    config_path: Path,
    fixture_path: Path,
    fixture_results: list[dict[str, Any]],
    volumes: CoarseVolumes,
    poro_inc_error: float,
    subvolume_metadata: list[dict[str, Any]],
) -> dict[str, Any]:
    fine_shape = tuple(int(value) for value in config["fine_shape_zyx"])
    spacing = [float(value) for value in config["coarse_spacing_xyz_m"]]
    if len(spacing) != 3 or any(value <= 0.0 for value in spacing):
        raise ValueError("coarse_spacing_xyz_m must contain three positive values in x, y, z order")
    inputs = {
        name: {"path": str(Path(config[name]).resolve()), "sha256": sha256(Path(config[name]))}
        for name in ("porosity_path", "facies_path", "poro_inc_path")
    }
    return {
        "schema_version": "stanford-vie-preprocessing.v2",
        "transform": "reshape(z,y,x)->mean_or_count_2x2x2->reverse_z_and_y",
        "axis_mapping": _axis_mapping(fine_shape, volumes.shape),
        "coarse_spacing_xyz_m": spacing,
        "strict_majority_rule": "channel_count>4 positive; channel_count<4 negative; channel_count==4 excluded",
        "channel_occupancy_definition": "channel_count / 8",
        "poro_inc_max_absolute_error": poro_inc_error,
        "inputs": inputs,
        "config": {"path": str(config_path.resolve()), "sha256": sha256(config_path)},
        "preprocessing_script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "co_location_fixture": {
            "path": str(fixture_path),
            "sha256": sha256(fixture_path),
            "results": fixture_results,
        },
        "global_counts": {
            "valid_voxels": int(np.count_nonzero(volumes.channel_label_valid)),
            "excluded_ties": int(np.count_nonzero(~volumes.channel_label_valid)),
        },
        "subvolumes": subvolume_metadata,
    }


def _validate_frozen_shapes(config: dict[str, Any], volumes: CoarseVolumes, subvolumes: list[dict[str, Any]]) -> None:
    expected_coarse = tuple(int(value) for value in config.get("expected_coarse_shape_zyx", volumes.shape))
    if volumes.shape != expected_coarse:
        raise RuntimeError(f"Coarse shape {volumes.shape} does not match the frozen shape {expected_coarse}")
    expected_subvolume = tuple(int(value) for value in config.get("expected_subvolume_shape_zyx", ()))
    if expected_subvolume and any(tuple(metadata["shape_zyx"]) != expected_subvolume for metadata in subvolumes):
        raise RuntimeError(f"A subvolume shape does not match the frozen shape {expected_subvolume}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create traceable co-located coarse Stanford VI-E porosity and channel products.")
    parser.add_argument("--config", type=Path, default=THIS_DIR / "preprocessing_config.json")
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output_root = Path(config["output_root"])
    output_root.mkdir(parents=True, exist_ok=True)
    volumes = coarse_volumes(config)
    _validate_frozen_shapes(config, volumes, [])
    fixtures, fixture_path = load_co_location_fixtures(config, config_path)
    fixture_results = validate_co_location_fixtures(
        volumes, fixtures, porosity_atol=float(config.get("fixture_porosity_absolute_tolerance", 1.0e-12))
    )
    reference = read_poro_inc(Path(config["poro_inc_path"]))
    if reference.size != volumes.porosity.size:
        raise RuntimeError("poro.inc value count does not match the coarse porosity grid")
    error = float(np.max(np.abs(volumes.porosity.reshape(-1) - reference)))
    if error > float(config["max_poro_inc_absolute_error"]):
        raise RuntimeError(f"Coarse porosity check failed: max absolute error {error:.9g}")
    subvolume_metadata: list[dict[str, Any]] = []
    for z_start, z_end in config["subvolumes"]:
        _, metadata = write_subvolume(output_root, volumes, int(z_start), int(z_end))
        subvolume_metadata.append(metadata)
    _validate_frozen_shapes(config, volumes, subvolume_metadata)
    metadata = build_manifest(
        config=config,
        config_path=config_path,
        fixture_path=fixture_path,
        fixture_results=fixture_results,
        volumes=volumes,
        poro_inc_error=error,
        subvolume_metadata=subvolume_metadata,
    )
    (output_root / "preprocessing_manifest.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
