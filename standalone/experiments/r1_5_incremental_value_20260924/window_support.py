from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from grid_io import GridData


@dataclass(frozen=True)
class WindowBatch:
    values: np.ndarray
    center_x: np.ndarray
    center_y: np.ndarray
    center_z: np.ndarray
    source_indices: np.ndarray

    @property
    def coordinates(self) -> pd.DataFrame:
        return pd.DataFrame({"x": self.center_x, "y": self.center_y, "z": self.center_z})


def extract_complete_windows(
    grid: GridData,
    window_size: tuple[int, int, int],
    move_step: tuple[int, int, int],
) -> WindowBatch:
    wx, wy, wz = (int(v) for v in window_size)
    sx, sy, sz = (int(v) for v in move_step)
    if wx <= 0 or wy <= 0 or wz not in (0, 1) and wz <= 0:
        raise ValueError("Window dimensions must be positive (2-D z may be 0 or 1)")
    if sx <= 0 or sy <= 0:
        raise ValueError("X and Y move steps must be positive")
    is_2d = grid.values.shape[0] == 1
    if is_2d:
        if wz not in (0, 1) or sz != 0:
            raise ValueError("2-D grids require window z=0/1 and step z=0")
        wz, sz = 1, 1
    elif wz <= 0 or sz <= 0:
        raise ValueError("3-D grids require positive z window and step")
    nz, ny, nx = grid.values.shape
    if wx > nx or wy > ny or wz > nz:
        raise ValueError("Full window does not fit inside the grid")

    records: list[np.ndarray] = []
    source_indices: list[np.ndarray] = []
    centers: list[tuple[int, int, int]] = []
    flat_indices = np.arange(grid.values.size).reshape(grid.values.shape)
    for z0 in range(0, nz - wz + 1, sz):
        for y0 in range(0, ny - wy + 1, sy):
            for x0 in range(0, nx - wx + 1, sx):
                support = grid.present[z0 : z0 + wz, y0 : y0 + wy, x0 : x0 + wx]
                if not support.all():
                    continue
                records.append(grid.values[z0 : z0 + wz, y0 : y0 + wy, x0 : x0 + wx].reshape(-1))
                source_indices.append(flat_indices[z0 : z0 + wz, y0 : y0 + wy, x0 : x0 + wx].reshape(-1))
                centers.append(
                    (
                        int(grid.x[x0 + (wx - 1) // 2]),
                        int(grid.y[y0 + (wy - 1) // 2]),
                        int(grid.z[z0 + (wz - 1) // 2]),
                    )
                )
    if not records:
        raise ValueError("No complete windows were extracted")
    center_array = np.asarray(centers, dtype=int)
    return WindowBatch(
        values=np.stack(records),
        center_x=center_array[:, 0],
        center_y=center_array[:, 1],
        center_z=center_array[:, 2],
        source_indices=np.stack(source_indices),
    )
