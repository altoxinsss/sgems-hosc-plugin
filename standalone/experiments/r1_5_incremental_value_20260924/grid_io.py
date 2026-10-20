from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class GridData:
    values: np.ndarray
    present: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    property_name: str
    source_path: Path
    title: str
    coordinates: pd.DataFrame

    @property
    def shape_xyz(self) -> tuple[int, int, int]:
        return (self.values.shape[2], self.values.shape[1], self.values.shape[0])


def _tokens(line: str) -> list[str]:
    return line.strip().split()


def _read_nonempty_lines(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Grid is too short: {path}")
    return lines


def _read_table(path: Path, property_name: str | None) -> tuple[str, pd.DataFrame, str]:
    lines = _read_nonempty_lines(path)
    title = lines[0].strip()
    second = _tokens(lines[1])
    if len(second) >= 2:
        # Legacy HOSC: ND NVAR, two comment lines, then index x y z value.
        nd = int(second[0])
        rows = []
        for line in lines[4 : 4 + nd]:
            fields = _tokens(line)
            if len(fields) < 5:
                raise ValueError(f"Malformed legacy row in {path}: {line!r}")
            rows.append((int(fields[1]), int(fields[2]), int(fields[3]), float(fields[4])))
        if len(rows) != nd:
            raise ValueError(f"Expected {nd} rows in {path}, found {len(rows)}")
        value_name = property_name or "value"
        return title, pd.DataFrame(rows, columns=["x", "y", "z", value_name]), value_name

    column_count = int(second[0])
    if column_count < 2 or len(lines) < 2 + column_count + 1:
        raise ValueError(f"Invalid columnar header in {path}")
    names = [lines[2 + index].strip() for index in range(column_count)]
    value_name = property_name or names[-1]
    lower_names = {name.lower(): name for name in names}
    if value_name.lower() not in lower_names:
        raise ValueError(f"Property {value_name!r} is not present in {path}; columns={names}")
    value_name = lower_names[value_name.lower()]
    rows = []
    for row_number, line in enumerate(lines[2 + column_count :], start=1):
        fields = _tokens(line)
        if not fields:
            continue
        if len(fields) < column_count:
            raise ValueError(f"Malformed row {row_number} in {path}: {line!r}")
        rows.append(fields[:column_count])
    frame = pd.DataFrame(rows, columns=names)
    coordinate_names = names[: names.index(value_name)]
    if len(coordinate_names) not in (1, 2, 3):
        raise ValueError(f"Could not infer 1-3 coordinate columns in {path}")
    normalized = pd.DataFrame()
    normalized["x"] = pd.to_numeric(frame[coordinate_names[0]], errors="raise").astype(int)
    normalized["y"] = (
        pd.to_numeric(frame[coordinate_names[1]], errors="raise").astype(int)
        if len(coordinate_names) >= 2
        else 0
    )
    normalized["z"] = (
        pd.to_numeric(frame[coordinate_names[2]], errors="raise").astype(int)
        if len(coordinate_names) >= 3
        else 0
    )
    normalized[value_name] = pd.to_numeric(frame[value_name], errors="raise").astype(float)
    return title, normalized, value_name


def read_grid(path: str | Path, property_name: str | None = None) -> GridData:
    source_path = Path(path).resolve()
    title, frame, value_name = _read_table(source_path, property_name)
    if frame.duplicated(["x", "y", "z"]).any():
        duplicates = frame.loc[frame.duplicated(["x", "y", "z"], keep=False), ["x", "y", "z"]]
        raise ValueError(f"Duplicate coordinates in {source_path}: {duplicates.head().to_dict('records')}")
    x = np.sort(frame["x"].unique())
    y = np.sort(frame["y"].unique())
    z = np.sort(frame["z"].unique())
    x_index = {int(value): index for index, value in enumerate(x)}
    y_index = {int(value): index for index, value in enumerate(y)}
    z_index = {int(value): index for index, value in enumerate(z)}
    values = np.full((len(z), len(y), len(x)), np.nan, dtype=float)
    present = np.zeros(values.shape, dtype=bool)
    for row in frame.itertuples(index=False):
        zi, yi, xi = z_index[int(row.z)], y_index[int(row.y)], x_index[int(row.x)]
        values[zi, yi, xi] = float(getattr(row, value_name))
        present[zi, yi, xi] = True
    coordinates = frame[["x", "y", "z"]].copy()
    coordinates[value_name] = frame[value_name].to_numpy(float)
    return GridData(values, present, x, y, z, value_name, source_path, title, coordinates)


def read_cli_output(path: str | Path) -> pd.DataFrame:
    grid = read_grid(path)
    frame = grid.coordinates.rename(columns={grid.property_name: "value"})
    return frame[["x", "y", "z", "value"]].sort_values(["z", "y", "x"]).reset_index(drop=True)


def write_map(
    path: str | Path,
    coordinates: pd.DataFrame | np.ndarray,
    values: Iterable[float],
    property_name: str,
    title: str = "R1-5 generated map",
) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    coords = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    data = np.asarray(list(values), dtype=float)
    if len(coords) != len(data):
        raise ValueError("Coordinate and value lengths differ")
    with target.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"{title}\n4\nx\ny\nz\n{property_name}\n")
        for row, value in zip(coords.itertuples(index=False), data):
            stream.write(f"{int(row.x)}\t{int(row.y)}\t{int(row.z)}\t{value:.17g}\n")


def align_maps(reference: pd.DataFrame, candidate: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    required = {"x", "y", "z", "value"}
    if not required.issubset(reference.columns) or not required.issubset(candidate.columns):
        raise ValueError("Maps must contain x, y, z, value columns")
    merged = reference.merge(candidate, on=["x", "y", "z"], how="outer", suffixes=("_reference", "_candidate"), indicator=True)
    if not (merged["_merge"] == "both").all():
        counts = merged["_merge"].value_counts().to_dict()
        raise ValueError(f"Coordinate supports differ: {counts}")
    merged = merged.sort_values(["z", "y", "x"])
    return merged["value_reference"].to_numpy(float), merged["value_candidate"].to_numpy(float)
