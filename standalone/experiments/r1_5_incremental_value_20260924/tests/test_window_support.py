from pathlib import Path

import numpy as np
import pandas as pd

from grid_io import read_grid, write_map
from window_support import extract_complete_windows


def test_complete_window_count_and_centers(tmp_path: Path):
    coordinates = pd.DataFrame([(x, y, 0) for y in range(1, 6) for x in range(1, 5)], columns=["x", "y", "z"])
    path = tmp_path / "grid.txt"
    write_map(path, coordinates, np.arange(20), "value")
    batch = extract_complete_windows(read_grid(path, "value"), (3, 3, 1), (1, 1, 0))
    assert batch.values.shape == (6, 9)
    assert list(batch.coordinates.iloc[0]) == [2, 2, 0]
    assert list(batch.coordinates.iloc[-1]) == [3, 4, 0]


def test_missing_cell_drops_affected_windows(tmp_path: Path):
    coordinates = pd.DataFrame([(x, y, 0) for y in range(4) for x in range(4) if (x, y) != (0, 0)], columns=["x", "y", "z"])
    path = tmp_path / "grid.txt"
    write_map(path, coordinates, np.arange(len(coordinates)), "value")
    batch = extract_complete_windows(read_grid(path, "value"), (2, 2, 1), (1, 1, 0))
    assert len(batch.values) == 8


def test_complete_window_count_includes_z_dimension(tmp_path: Path):
    coordinates = pd.DataFrame(
        [(x, y, z) for z in range(4) for y in range(4) for x in range(5)], columns=["x", "y", "z"]
    )
    path = tmp_path / "volume.txt"
    write_map(path, coordinates, np.arange(len(coordinates)), "value")
    batch = extract_complete_windows(read_grid(path, "value"), (3, 2, 2), (1, 2, 1))
    assert batch.values.shape == (18, 12)
    assert list(batch.coordinates.iloc[0]) == [1, 0, 0]
    assert list(batch.coordinates.iloc[-1]) == [3, 2, 2]
