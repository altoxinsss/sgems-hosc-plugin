from __future__ import annotations

from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

from grid_io import GridData


def read_facies_dat(
    path: str | Path,
    shape_zyx: tuple[int, int, int],
    header_lines: int = 3,
) -> np.ndarray:
    """Read the Stanford VI-E integer facies volume in z/y/x order."""
    target = Path(path)
    expected = int(np.prod(shape_zyx))
    with target.open("r", encoding="utf-8") as stream:
        for _ in range(header_lines):
            if not stream.readline():
                raise ValueError(f"Facies file ends before its {header_lines} header lines: {target}")
        values = np.fromfile(stream, dtype=np.int16, sep=" ")
    if len(values) != expected:
        raise ValueError(f"Facies value count {len(values)} does not match shape {shape_zyx} ({expected})")
    return values.reshape(shape_zyx)


def downsample_facies_mode(facies: np.ndarray, target_shape_zyx: tuple[int, int, int]) -> np.ndarray:
    """Downsample categorical facies by blockwise mode, never numeric averaging."""
    fine = np.asarray(facies)
    if fine.ndim != 3:
        raise ValueError("facies must be a three-dimensional z/y/x array")
    target_shape = tuple(int(value) for value in target_shape_zyx)
    if any(value <= 0 for value in target_shape) or any(fine.shape[index] % target_shape[index] for index in range(3)):
        raise ValueError(f"Facies shape {fine.shape} is not an integer block multiple of {target_shape}")
    factors = tuple(fine.shape[index] // target_shape[index] for index in range(3))
    blocks = fine.reshape(
        target_shape[0], factors[0], target_shape[1], factors[1], target_shape[2], factors[2]
    )
    categories = np.unique(fine)
    counts = np.stack([np.sum(blocks == category, axis=(1, 3, 5)) for category in categories])
    return categories[np.argmax(counts, axis=0)]


def block_reduce_mean(values: np.ndarray, target_shape_zyx: tuple[int, int, int]) -> np.ndarray:
    """Block-average a continuous z/y/x volume using an exact integer ratio."""
    fine = np.asarray(values, dtype=float)
    target = tuple(int(value) for value in target_shape_zyx)
    if fine.ndim != 3 or any(value <= 0 for value in target):
        raise ValueError("values and target_shape_zyx must describe positive 3-D arrays")
    if any(fine.shape[axis] % target[axis] for axis in range(3)):
        raise ValueError(f"Shape {fine.shape} is not an integer block multiple of {target}")
    factors = tuple(fine.shape[axis] // target[axis] for axis in range(3))
    return fine.reshape(target[0], factors[0], target[1], factors[1], target[2], factors[2]).mean(axis=(1, 3, 5))


def channel_count_2x2x2(facies: np.ndarray, channel_code: int = 2) -> np.ndarray:
    """Return the channel occupancy count in every exact 2x2x2 fine block."""
    fine = np.asarray(facies)
    if fine.ndim != 3 or any(length % 2 for length in fine.shape):
        raise ValueError("facies must have even z/y/x dimensions for 2x2x2 occupancy")
    target = tuple(length // 2 for length in fine.shape)
    blocks = fine.reshape(target[0], 2, target[1], 2, target[2], 2)
    return np.sum(blocks == int(channel_code), axis=(1, 3, 5), dtype=np.int8)


def strict_majority_channel_labels(channel_count: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return binary labels and a validity mask; 4/8 ties are excluded."""
    counts = np.asarray(channel_count, dtype=np.int8)
    if np.any((counts < 0) | (counts > 8)):
        raise ValueError("channel_count must lie in [0, 8]")
    valid = counts != 4
    labels = counts > 4
    return labels, valid


def channel_facies_mask(facies: np.ndarray, channel_code: int = 2) -> np.ndarray:
    return np.asarray(facies) == int(channel_code)


def subset_facies_by_grid_coordinates(
    facies: np.ndarray,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    z_coordinates: np.ndarray,
) -> np.ndarray:
    """Select categorical facies with the one-based Stanford grid coordinates."""
    coarse = np.asarray(facies)
    if coarse.ndim != 3:
        raise ValueError("facies must be a three-dimensional z/y/x array")
    x_index = np.asarray(x_coordinates, dtype=int) - 1
    y_index = np.asarray(y_coordinates, dtype=int) - 1
    z_index = np.asarray(z_coordinates, dtype=int) - 1
    for name, indices, size in (("x", x_index, coarse.shape[2]), ("y", y_index, coarse.shape[1]), ("z", z_index, coarse.shape[0])):
        if len(indices) == 0 or np.any(indices < 0) or np.any(indices >= size):
            raise ValueError(f"Grid {name} coordinates are outside the facies domain 1..{size}")
    return coarse[np.ix_(z_index, y_index, x_index)]


def make_high_porosity_mask(grid: GridData, threshold: float) -> np.ndarray:
    return grid.present & (grid.values >= float(threshold))


def _neighbors(z: int, y: int, x: int, shape: tuple[int, int, int]):
    nz, ny, nx = shape
    if x > 0:
        yield z, y, x - 1
    if x + 1 < nx:
        yield z, y, x + 1
    if y > 0:
        yield z, y - 1, x
    if y + 1 < ny:
        yield z, y + 1, x
    if nz > 1 and z > 0:
        yield z - 1, y, x
    if nz > 1 and z + 1 < nz:
        yield z + 1, y, x


def make_main_channel_mask(mask: np.ndarray, component_count: int) -> np.ndarray:
    source = np.asarray(mask, dtype=bool)
    if component_count <= 0:
        raise ValueError("component_count must be positive")
    visited = np.zeros_like(source, dtype=bool)
    components: list[list[tuple[int, int, int]]] = []
    for start in zip(*np.where(source & ~visited)):
        if visited[start]:
            continue
        queue = deque([start])
        visited[start] = True
        component: list[tuple[int, int, int]] = []
        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbor in _neighbors(*current, source.shape):
                if source[neighbor] and not visited[neighbor]:
                    visited[neighbor] = True
                    queue.append(neighbor)
        components.append(component)
    components.sort(key=lambda item: (-len(item), min(item) if item else (0, 0, 0)))
    result = np.zeros_like(source, dtype=bool)
    for component in components[:component_count]:
        for index in component:
            result[index] = True
    return result


def labels_at_window_centers(mask: np.ndarray, centers: pd.DataFrame, grid: GridData) -> np.ndarray:
    x_index = {int(value): index for index, value in enumerate(grid.x)}
    y_index = {int(value): index for index, value in enumerate(grid.y)}
    z_index = {int(value): index for index, value in enumerate(grid.z)}
    labels = []
    for row in centers.itertuples(index=False):
        try:
            labels.append(bool(mask[z_index[int(row.z)], y_index[int(row.y)], x_index[int(row.x)]]))
        except KeyError as error:
            raise ValueError(f"Window center is outside the source grid: {row}") from error
    return np.asarray(labels, dtype=int)
