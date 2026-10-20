from pathlib import Path

import numpy as np
import pandas as pd

from grid_io import align_maps, read_grid, write_map


def test_columnar_round_trip(tmp_path: Path):
    path = tmp_path / "grid.txt"
    coordinates = pd.DataFrame({"x": [1, 2, 1, 2], "y": [1, 1, 2, 2], "z": [0, 0, 0, 0]})
    values = np.array([0.1, 0.2, 0.3, 0.4])
    write_map(path, coordinates, values, "value")
    grid = read_grid(path, "value")
    assert grid.shape_xyz == (2, 2, 1)
    assert grid.present.all()
    assert np.allclose(grid.values.ravel(), values)
    frame = coordinates.copy()
    frame["value"] = values
    reference, candidate = align_maps(frame, frame.sample(frac=1.0, random_state=1))
    assert np.array_equal(reference, candidate)


def test_duplicate_coordinates_are_rejected(tmp_path: Path):
    path = tmp_path / "duplicate.grid"
    path.write_text("g\n3\nx\ny\nv\n1 1 0.1\n1 1 0.2\n", encoding="utf-8")
    try:
        read_grid(path, "v")
    except ValueError as error:
        assert "Duplicate" in str(error)
    else:
        raise AssertionError("duplicate coordinates were accepted")
