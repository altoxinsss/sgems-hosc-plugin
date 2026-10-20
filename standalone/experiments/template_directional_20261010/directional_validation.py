from __future__ import annotations

"""Numerical utilities used by the independent Section 4.4 experiment.

The production ``template_cli`` accepts directional lag descriptions, whereas
the shape checks in this module use explicit offsets. Both paths are kept
separate in the output manifest so a reference calculation is never reported
as a production CLI result.
"""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class DenseGrid:
    """A dense grid stored in z, y, x order with physical coordinate vectors."""

    values: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    value_name: str

    @property
    def is_2d(self) -> bool:
        return self.values.shape[0] == 1


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rotate_offsets_90(offsets: np.ndarray) -> np.ndarray:
    offsets = np.asarray(offsets, dtype=int)
    if offsets.ndim != 2 or offsets.shape[1] != 2:
        raise ValueError("offsets must have shape (points, 2)")
    return np.column_stack((-offsets[:, 1], offsets[:, 0]))


def directional_offset(angle_degrees: int, lag: int) -> tuple[int, int]:
    directions = {0: (1, 0), 45: (1, 1), 90: (0, 1), 135: (-1, 1)}
    if angle_degrees not in directions or lag <= 0:
        raise ValueError("angle must be 0, 45, 90, or 135 and lag must be positive")
    dx, dy = directions[angle_degrees]
    return dx * lag, dy * lag


def transform_offsets_2d(offsets: np.ndarray, angle_degrees: int, lag: int) -> np.ndarray:
    """Rotate a template expressed in the x/y basis onto an integer direction."""

    if lag <= 0:
        raise ValueError("lag must be positive")
    bases = {
        0: np.array(((1, 0), (0, 1)), dtype=int),
        # Each matrix maps a local (along-direction, transverse-direction)
        # row vector to an x/y grid offset through ``local @ matrix.T``.
        45: np.array(((1, -1), (1, 1)), dtype=int),
        90: np.array(((0, -1), (1, 0)), dtype=int),
        135: np.array(((-1, -1), (1, -1)), dtype=int),
    }
    if angle_degrees not in bases:
        raise ValueError("angle must be 0, 45, 90, or 135")
    source = np.asarray(offsets, dtype=int)
    if source.ndim != 2 or source.shape[1] != 2:
        raise ValueError("offsets must have shape (points, 2)")
    return (source @ bases[angle_degrees].T) * lag


def fourth_joint_cumulant(values: np.ndarray) -> float:
    """Evaluate the joint fourth cumulant for exactly four variables."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("values must have shape (samples, 4)")
    if array.shape[0] == 0:
        return float("nan")
    centered = array - array.mean(axis=0, keepdims=True)
    m4 = np.mean(np.prod(centered, axis=1))
    pairings = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))
    correction = sum(
        np.mean(centered[:, a] * centered[:, b]) * np.mean(centered[:, c] * centered[:, d])
        for a, b, c, d in pairings
    )
    return float(m4 - correction)


def template_response(field: np.ndarray, offsets: np.ndarray) -> float:
    """Return a fourth-order response for a 2D template over valid translations."""

    image = np.asarray(field, dtype=float)
    if image.ndim != 2:
        raise ValueError("field must be two-dimensional")
    xy = np.asarray(offsets, dtype=int)
    if xy.shape != (4, 2):
        raise ValueError("a fourth-order 2D template must contain four x/y offsets")
    xyz = np.column_stack((xy, np.zeros(4, dtype=int)))
    return template_response_nd(image[np.newaxis, :, :], xyz)


def template_response_nd(volume: np.ndarray, offsets_xyz: np.ndarray) -> float:
    """Return a fourth-order response for x/y/z offsets over valid translations."""

    values = np.asarray(volume, dtype=float)
    offsets = np.asarray(offsets_xyz, dtype=int)
    if values.ndim != 3:
        raise ValueError("volume must use z, y, x ordering")
    if offsets.shape != (4, 3):
        raise ValueError("a fourth-order 3D template must contain four x/y/z offsets")

    minimum = offsets.min(axis=0)
    maximum = offsets.max(axis=0)
    nx, ny, nz = values.shape[2], values.shape[1], values.shape[0]
    if np.any(maximum - minimum + 1 > np.array((nx, ny, nz))):
        return float("nan")

    x_start, x_stop = -minimum[0], nx - maximum[0]
    y_start, y_stop = -minimum[1], ny - maximum[1]
    z_start, z_stop = -minimum[2], nz - maximum[2]
    translated = []
    for dx, dy, dz in offsets:
        translated.append(
            values[z_start + dz : z_stop + dz, y_start + dy : y_stop + dy, x_start + dx : x_stop + dx].reshape(-1)
        )
    sample = np.column_stack(translated)
    return fourth_joint_cumulant(sample[np.isfinite(sample).all(axis=1)])


def iaaft_surrogate(field: np.ndarray, seed: int, iterations: int = 200) -> np.ndarray:
    """Generate a histogram-preserving Fourier-amplitude surrogate without SciPy."""

    source = np.asarray(field, dtype=float)
    if source.ndim not in (2, 3) or not np.isfinite(source).all():
        raise ValueError("IAAFT requires a finite 2D or 3D field")
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    rng = np.random.default_rng(seed)
    sorted_values = np.sort(source.reshape(-1))
    target_amplitude = np.abs(np.fft.fftn(source))
    current = rng.permutation(source.reshape(-1)).reshape(source.shape)
    for _ in range(iterations):
        phase = np.angle(np.fft.fftn(current))
        constrained = np.fft.ifftn(target_amplitude * np.exp(1j * phase)).real
        ranks = np.argsort(np.argsort(constrained.reshape(-1), kind="mergesort"), kind="mergesort")
        current = sorted_values[ranks].reshape(source.shape)
    return current


def empirical_variogram(field: np.ndarray, lags: Sequence[int]) -> np.ndarray:
    """Compute the axis-averaged semivariogram for a finite 2D or 3D field."""

    array = np.asarray(field, dtype=float)
    if array.ndim not in (2, 3) or not np.isfinite(array).all():
        raise ValueError("variogram requires a finite 2D or 3D field")
    result: list[float] = []
    for lag in lags:
        if lag <= 0:
            raise ValueError("variogram lags must be positive")
        samples = []
        for axis, length in enumerate(array.shape):
            if lag >= length:
                continue
            lower = [slice(None)] * array.ndim
            upper = [slice(None)] * array.ndim
            lower[axis], upper[axis] = slice(None, -lag), slice(lag, None)
            difference = array[tuple(upper)] - array[tuple(lower)]
            samples.append(np.mean(difference * difference) / 2.0)
        result.append(float(np.mean(samples)) if samples else float("nan"))
    return np.asarray(result, dtype=float)


def variogram_relative_rmse(reference: np.ndarray, candidate: np.ndarray) -> float:
    reference, candidate = np.asarray(reference, dtype=float), np.asarray(candidate, dtype=float)
    finite = np.isfinite(reference) & np.isfinite(candidate)
    if not finite.any():
        return float("inf")
    scale = max(float(np.sqrt(np.mean(reference[finite] ** 2))), np.finfo(float).eps)
    return float(np.sqrt(np.mean((reference[finite] - candidate[finite]) ** 2)) / scale)


def validate_iaaft_surrogate(reference: np.ndarray, candidate: np.ndarray, lags: Sequence[int],
                             max_variogram_relative_rmse: float) -> dict[str, float | bool]:
    """Return the histogram and variogram gates used for every surrogate."""

    same_histogram = bool(np.array_equal(np.sort(reference.reshape(-1)), np.sort(candidate.reshape(-1))))
    relative_rmse = variogram_relative_rmse(empirical_variogram(reference, lags), empirical_variogram(candidate, lags))
    return {
        "histogram_exact": same_histogram,
        "variogram_relative_rmse": relative_rmse,
        "variogram_pass": bool(relative_rmse <= max_variogram_relative_rmse),
        "accepted": bool(same_histogram and relative_rmse <= max_variogram_relative_rmse),
    }


def generate_accepted_iaaft_surrogates(field: np.ndarray, count: int, seed: int, iterations: int,
                                       lags: Sequence[int], max_variogram_relative_rmse: float,
                                       max_attempts: int) -> tuple[list[np.ndarray], list[dict[str, float | bool | int]]]:
    """Generate an exact number of accepted IAAFT fields or fail explicitly."""

    if count <= 0 or max_attempts < count:
        raise ValueError("count must be positive and max_attempts must be at least count")
    accepted: list[np.ndarray] = []
    audit: list[dict[str, float | bool | int]] = []
    for attempt in range(max_attempts):
        candidate = iaaft_surrogate(field, seed + attempt, iterations)
        verdict = validate_iaaft_surrogate(field, candidate, lags, max_variogram_relative_rmse)
        record: dict[str, float | bool | int] = {"attempt": attempt + 1, "seed": seed + attempt, **verdict}
        audit.append(record)
        if bool(verdict["accepted"]):
            accepted.append(candidate)
        if len(accepted) == count:
            return accepted, audit
    raise RuntimeError(f"IAAFT acceptance gate produced {len(accepted)} accepted fields after {max_attempts} attempts; expected {count}.")


def _component_sizes(mask: np.ndarray) -> list[int]:
    binary = np.asarray(mask, dtype=bool)
    visited = np.zeros(binary.shape, dtype=bool)
    neighbor_offsets = []
    for axis in range(binary.ndim):
        for delta in (-1, 1):
            vector = [0] * binary.ndim
            vector[axis] = delta
            neighbor_offsets.append(tuple(vector))
    sizes: list[int] = []
    for origin in np.argwhere(binary):
        start = tuple(int(value) for value in origin)
        if visited[start]:
            continue
        stack, size = [start], 0
        visited[start] = True
        while stack:
            current = stack.pop()
            size += 1
            for offset in neighbor_offsets:
                neighbor = tuple(current[index] + offset[index] for index in range(binary.ndim))
                if any(value < 0 or value >= binary.shape[index] for index, value in enumerate(neighbor)):
                    continue
                if binary[neighbor] and not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(neighbor)
        sizes.append(size)
    return sizes


def directional_run_lengths(mask: np.ndarray, axis_labels: Sequence[str] | None = None) -> dict[str, float | int]:
    """Summarize occupied-run lengths along every grid axis."""

    binary = np.asarray(mask, dtype=bool)
    labels = tuple(axis_labels) if axis_labels is not None else tuple(f"axis_{axis}" for axis in range(binary.ndim))
    if len(labels) != binary.ndim:
        raise ValueError("axis_labels must have one label per array axis")
    result: dict[str, float | int] = {}
    for axis, label in enumerate(labels):
        lines = np.moveaxis(binary, axis, -1).reshape(-1, binary.shape[axis])
        lengths: list[int] = []
        for line in lines:
            padded = np.concatenate(([False], line, [False])).astype(int)
            transitions = np.flatnonzero(np.diff(padded))
            lengths.extend((transitions[1::2] - transitions[::2]).tolist())
        result[f"run_count_{label}"] = len(lengths)
        result[f"run_max_{label}"] = int(max(lengths)) if lengths else 0
        result[f"run_mean_{label}"] = float(np.mean(lengths)) if lengths else 0.0
    return result


def connectivity_metrics(mask: np.ndarray, axis_labels: Sequence[str] | None = None) -> dict[str, float | int]:
    """Return face-connected component and directional-run diagnostics."""

    binary = np.asarray(mask, dtype=bool)
    if binary.ndim not in (2, 3):
        raise ValueError("connectivity diagnostics support 2D or 3D masks")
    active = int(binary.sum())
    components = _component_sizes(binary)
    largest = max(components) if components else 0
    result: dict[str, float | int] = {
        "total_voxels": int(binary.size), "active_voxels": active,
        "active_fraction": float(active / binary.size), "component_count": len(components),
        "largest_component_voxels": int(largest),
        "largest_component_fraction_of_active": float(largest / active) if active else 0.0,
    }
    labels = axis_labels if axis_labels is not None else (("y", "x") if binary.ndim == 2 else ("z", "y", "x"))
    result.update(directional_run_lengths(binary, labels))
    return result


def read_columnar_grid(path: Path, property_name: str | None = None) -> DenseGrid:
    """Read the named-column grid format emitted by the production CLI."""

    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 6:
        raise ValueError(f"grid is too short: {path}")
    try:
        column_count = int(lines[1].strip())
    except ValueError as error:
        raise ValueError(f"only named-column grids are supported: {path}") from error
    headers = [line.strip() for line in lines[2 : 2 + column_count]]
    data_start = 2 + column_count
    records = np.loadtxt(path, dtype=float, skiprows=data_start, ndmin=2)
    if records.shape[1] != column_count:
        raise ValueError(f"grid data column count is invalid: {path}")
    lowered = [header.lower() for header in headers]
    try:
        x_index, y_index = lowered.index("x"), lowered.index("y")
    except ValueError as error:
        raise ValueError(f"grid must contain x and y coordinates: {path}") from error
    z_index = lowered.index("z") if "z" in lowered else None
    candidates = [index for index in range(column_count) if index not in {x_index, y_index, z_index}]
    if property_name is None:
        if len(candidates) != 1:
            raise ValueError(f"grid property is ambiguous: {path}")
        value_index = candidates[0]
    else:
        try:
            value_index = lowered.index(property_name.lower())
        except ValueError as error:
            raise ValueError(f"property '{property_name}' is not present in {path}") from error
    x, y = np.unique(records[:, x_index].astype(int)), np.unique(records[:, y_index].astype(int))
    z = np.unique(records[:, z_index].astype(int)) if z_index is not None else np.array([0], dtype=int)
    output = np.full((z.size, y.size, x.size), np.nan, dtype=float)
    x_lookup, y_lookup, z_lookup = ({int(value): index for index, value in enumerate(axis)} for axis in (x, y, z))
    for row in records:
        iz = z_lookup[int(row[z_index])] if z_index is not None else 0
        output[iz, y_lookup[int(row[y_index])], x_lookup[int(row[x_index])]] = row[value_index]
    if not np.isfinite(output).all():
        raise ValueError(f"grid is not a complete rectilinear volume: {path}")
    return DenseGrid(output, x, y, z, headers[value_index])


def write_columnar_grid(path: Path, values: np.ndarray, value_name: str, *, x: Iterable[int] | None = None,
                        y: Iterable[int] | None = None, z: Iterable[int] | None = None,
                        allow_nan: bool = False) -> None:
    """Write a complete finite 2D or 3D named-column grid for ``template_cli``."""

    array = np.asarray(values, dtype=float)
    if array.ndim not in (2, 3) or (not allow_nan and not np.isfinite(array).all()):
        raise ValueError("a finite 2D or 3D array is required unless allow_nan is enabled")
    if allow_nan and np.isinf(array).any():
        raise ValueError("allow_nan permits missing values but not infinities")
    path.parent.mkdir(parents=True, exist_ok=True)
    if array.ndim == 2:
        ny, nx = array.shape
        xx, yy = np.asarray(list(x) if x is not None else np.arange(nx), dtype=int), np.asarray(list(y) if y is not None else np.arange(ny), dtype=int)
        if xx.size != nx or yy.size != ny:
            raise ValueError("2D coordinate dimensions do not match the data")
        rows, header = ((int(xx[ix]), int(yy[iy]), float(array[iy, ix])) for iy in range(ny) for ix in range(nx)), ("x", "y", value_name)
    else:
        nz, ny, nx = array.shape
        xx, yy, zz = (np.asarray(list(axis) if axis is not None else np.arange(length), dtype=int) for axis, length in ((x, nx), (y, ny), (z, nz)))
        if xx.size != nx or yy.size != ny or zz.size != nz:
            raise ValueError("3D coordinate dimensions do not match the data")
        rows = ((int(xx[ix]), int(yy[iy]), int(zz[iz]), float(array[iz, iy, ix])) for iz in range(nz) for iy in range(ny) for ix in range(nx))
        header = ("x", "y", "z", value_name)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"{value_name}\n{len(header)}\n")
        stream.write("\n".join(header) + "\n")
        for row in rows:
            stream.write("\t".join(format(value, ".17g") for value in row) + "\n")


def rotation_exchange_error(responses: dict[int, float], rotated_responses: dict[int, float]) -> dict[str, float]:
    """Measure the expected 0<->90 and 45<->135 response exchange."""

    pairs = ((0, 90), (45, 135))
    return {f"{source}_to_{target}": abs(responses[source] - rotated_responses[target]) for source, target in pairs}


# Formal Section 4.4 validation -------------------------------------------------
#
# The helpers below deliberately know nothing about cumulant or connectivity
# values.  Surrogate acceptance is based only on the frozen H2 offset set and
# the field's histogram, preventing post-selection on a higher-order result.

_SIGNED_2D_DIRECTIONS: dict[int, tuple[int, int]] = {
    0: (1, 0),
    45: (1, 1),
    90: (0, 1),
    135: (-1, 1),
    180: (-1, 0),
    225: (-1, -1),
    270: (0, -1),
    315: (1, -1),
}


def signed_integer_direction(angle_degrees: int) -> tuple[int, int]:
    """Return a signed grid direction without reducing angles modulo 180."""

    try:
        return _SIGNED_2D_DIRECTIONS[int(angle_degrees)]
    except KeyError as error:
        raise ValueError("angle must be one of 0,45,...,315 degrees") from error


def _as_offset_array(offsets: Sequence[Sequence[int]], dimensions: int) -> np.ndarray:
    array = np.asarray(offsets, dtype=int)
    if array.ndim != 2 or array.shape[1] != dimensions:
        raise ValueError(f"offsets must have shape (count, {dimensions})")
    return array


def semivariogram_at_offsets(
    field: np.ndarray,
    offsets: Sequence[Sequence[int]],
    valid_mask: np.ndarray | None = None,
) -> list[dict[str, float | int | tuple[int, ...]]]:
    """Calculate masked semivariograms and pair counts for signed offsets."""

    values = np.asarray(field, dtype=float)
    if values.ndim not in (2, 3):
        raise ValueError("field must be two- or three-dimensional")
    vector_array = _as_offset_array(offsets, values.ndim)
    mask = np.isfinite(values) if valid_mask is None else np.asarray(valid_mask, dtype=bool) & np.isfinite(values)
    if mask.shape != values.shape:
        raise ValueError("valid_mask shape must equal field shape")

    results: list[dict[str, float | int | tuple[int, ...]]] = []
    for vector in vector_array:
        if not np.any(vector):
            valid = mask
            difference = np.zeros_like(values, dtype=float)
        else:
            lower_slices: list[slice] = []
            upper_slices: list[slice] = []
            for axis, shift in enumerate(vector.tolist()):
                if abs(shift) >= values.shape[axis]:
                    lower_slices = []
                    break
                lower_slices.append(slice(max(0, -shift), values.shape[axis] - max(0, shift)))
                upper_slices.append(slice(max(0, shift), values.shape[axis] - max(0, -shift)))
            if not lower_slices:
                results.append({"offset": tuple(int(item) for item in vector), "gamma": float("nan"), "pair_count": 0})
                continue
            low = tuple(lower_slices)
            high = tuple(upper_slices)
            valid = mask[low] & mask[high]
            difference = values[high] - values[low]
        pair_count = int(np.count_nonzero(valid))
        gamma = float(np.mean(np.square(difference[valid])) / 2.0) if pair_count else float("nan")
        results.append({"offset": tuple(int(item) for item in vector), "gamma": gamma, "pair_count": pair_count})
    return results


def directional_semivariogram_2d(
    field: np.ndarray,
    lags: Sequence[int],
    valid_mask: np.ndarray | None = None,
    angles: Sequence[int] = (0, 45, 90, 135),
) -> list[dict[str, float | int | tuple[int, ...]]]:
    """Return one row for each signed 2D direction and positive lag."""

    image = np.asarray(field, dtype=float)
    if image.ndim != 2:
        raise ValueError("directional_semivariogram_2d requires a 2D field")
    rows: list[dict[str, float | int | tuple[int, ...]]] = []
    for angle in angles:
        dx, dy = signed_integer_direction(int(angle))
        for lag in lags:
            if int(lag) <= 0:
                raise ValueError("directional variogram lags must be positive")
            record = semivariogram_at_offsets(image, [(dy * int(lag), dx * int(lag))], valid_mask)[0]
            rows.append({"angle_degrees": int(angle), "lag_cells": int(lag), **record})
    return rows


def axis_semivariogram_3d(
    volume: np.ndarray,
    lags: Sequence[int],
    valid_mask: np.ndarray | None = None,
) -> list[dict[str, float | int | tuple[int, ...]]]:
    """Return semivariograms along the physical X/Y/Z grid axes."""

    values = np.asarray(volume, dtype=float)
    if values.ndim != 3:
        raise ValueError("axis_semivariogram_3d requires a 3D field")
    directions = {"x": (0, 0, 1), "y": (0, 1, 0), "z": (1, 0, 0)}
    rows: list[dict[str, float | int | tuple[int, ...]]] = []
    for name, vector in directions.items():
        for lag in lags:
            if int(lag) <= 0:
                raise ValueError("axis variogram lags must be positive")
            offset = tuple(int(lag) * item for item in vector)
            record = semivariogram_at_offsets(values, [offset], valid_mask)[0]
            rows.append({"axis": name, "lag_cells": int(lag), **record})
    return rows


def h2_offsets_from_templates(templates: Iterable[np.ndarray]) -> np.ndarray:
    """Freeze all signed pair offsets implied by a collection of templates."""

    values: set[tuple[int, ...]] = set()
    dimensions: int | None = None
    for template in templates:
        points = np.asarray(template, dtype=int)
        if points.ndim != 2 or points.shape[0] < 2:
            raise ValueError("every template must contain at least two points")
        if dimensions is None:
            dimensions = int(points.shape[1])
        if points.shape[1] != dimensions:
            raise ValueError("all templates must have the same dimensionality")
        for left in points:
            for right in points:
                values.add(tuple(int(item) for item in right - left))
    if dimensions is None:
        raise ValueError("at least one template is required")
    # Zero lag is useful for an audit but does not define a directional pair
    # match.  It remains in H2 output and receives a deterministic gamma=0.
    return np.asarray(sorted(values), dtype=int).reshape(-1, dimensions)


def global_variance(field: np.ndarray, valid_mask: np.ndarray | None = None) -> float:
    values = np.asarray(field, dtype=float)
    mask = np.isfinite(values) if valid_mask is None else np.asarray(valid_mask, dtype=bool) & np.isfinite(values)
    if mask.shape != values.shape:
        raise ValueError("valid_mask shape must equal field shape")
    if int(mask.sum()) < 2:
        return float("nan")
    return float(np.var(values[mask], ddof=0))


def h2_match_diagnostics(
    reference: np.ndarray,
    candidate: np.ndarray,
    h2_offsets: Sequence[Sequence[int]],
    valid_mask: np.ndarray | None = None,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> dict[str, object]:
    """Compare all frozen H2 offsets using the stated variance-normalized gate."""

    ref = np.asarray(reference, dtype=float)
    other = np.asarray(candidate, dtype=float)
    if ref.shape != other.shape:
        raise ValueError("reference and candidate shapes must match")
    variance = global_variance(ref, valid_mask)
    if not np.isfinite(variance) or variance <= 0.0:
        return {"accepted": False, "reason": "reference_variance_nonpositive", "rows": []}
    reference_rows = semivariogram_at_offsets(ref, h2_offsets, valid_mask)
    candidate_rows = semivariogram_at_offsets(other, h2_offsets, valid_mask)
    rows: list[dict[str, object]] = []
    weighted_error = 0.0
    total_pairs = 0
    maximum = 0.0
    valid = True
    for baseline, actual in zip(reference_rows, candidate_rows, strict=True):
        pair_count = min(int(baseline["pair_count"]), int(actual["pair_count"]))
        delta = abs(float(actual["gamma"]) - float(baseline["gamma"]))
        relative = delta / variance if pair_count and np.isfinite(delta) else float("inf")
        maximum = max(maximum, relative)
        weighted_error += pair_count * relative * relative
        total_pairs += pair_count
        valid = valid and pair_count > 0 and np.isfinite(relative)
        rows.append({
            "offset": baseline["offset"],
            "reference_gamma": baseline["gamma"],
            "candidate_gamma": actual["gamma"],
            "pair_count": pair_count,
            "absolute_relative_error": relative,
        })
    rmse = float(np.sqrt(weighted_error / total_pairs)) if total_pairs else float("inf")
    accepted = bool(valid and maximum <= max_absolute_relative_error and rmse <= max_weighted_relative_rmse)
    return {
        "accepted": accepted,
        "rows": rows,
        "variance": variance,
        "maximum_absolute_relative_error": maximum,
        "weighted_relative_rmse": rmse,
        "pair_count_total": total_pairs,
    }


def _swap_objective(
    reference: np.ndarray,
    candidate: np.ndarray,
    h2_offsets: np.ndarray,
    valid_mask: np.ndarray | None,
    *,
    max_absolute_relative_error: float,
    max_weighted_relative_rmse: float,
) -> tuple[float, tuple[float, float, float]]:
    diagnostic = h2_match_diagnostics(reference, candidate, h2_offsets, valid_mask,
                                      max_absolute_relative_error=float("inf"),
                                      max_weighted_relative_rmse=float("inf"))
    rmse = float(diagnostic["weighted_relative_rmse"])
    maximum = float(diagnostic["maximum_absolute_relative_error"])
    if not (np.isfinite(rmse) and np.isfinite(maximum)):
        return float("inf"), (float("inf"), float("inf"), float("inf"))
    # Rank candidates by frozen-gate violation before optimizing within the
    # feasible set.  The scalar score remains only for the annealing proposal;
    # the lexicographic rank decides the result returned to the caller.
    violation = max(0.0, maximum - max_absolute_relative_error) + max(0.0, rmse - max_weighted_relative_rmse)
    return rmse, (violation, rmse, maximum)


def optimize_h2_by_swaps(
    reference: np.ndarray,
    candidate: np.ndarray,
    h2_offsets: Sequence[Sequence[int]],
    seed: int,
    iterations: int,
    valid_mask: np.ndarray | None = None,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> np.ndarray:
    """Use deterministic simulated swaps while preserving the value multiset."""

    if iterations <= 0:
        return np.asarray(candidate, dtype=float).copy()
    target = np.asarray(reference, dtype=float)
    current = np.asarray(candidate, dtype=float).copy()
    mask = np.isfinite(current) if valid_mask is None else np.asarray(valid_mask, dtype=bool) & np.isfinite(current)
    positions = np.flatnonzero(mask.reshape(-1))
    if positions.size < 2:
        return current
    offsets = _as_offset_array(h2_offsets, current.ndim)
    rng = np.random.default_rng(seed)
    score, quality = _swap_objective(
        target, current, offsets, mask,
        max_absolute_relative_error=max_absolute_relative_error,
        max_weighted_relative_rmse=max_weighted_relative_rmse,
    )
    best = current.copy()
    best_quality = quality
    for iteration in range(iterations):
        first, second = rng.choice(positions, size=2, replace=False)
        flat = current.reshape(-1)
        if flat[first] == flat[second]:
            continue
        flat[first], flat[second] = flat[second], flat[first]
        candidate_score, candidate_quality = _swap_objective(
            target, current, offsets, mask,
            max_absolute_relative_error=max_absolute_relative_error,
            max_weighted_relative_rmse=max_weighted_relative_rmse,
        )
        temperature = max(1.0e-9, 0.05 * (1.0 - iteration / max(1, iterations)))
        accept = candidate_score <= score or rng.random() < np.exp((score - candidate_score) / temperature)
        if accept:
            score = candidate_score
            if candidate_quality < best_quality:
                best = current.copy()
                best_quality = candidate_quality
        else:
            flat[first], flat[second] = flat[second], flat[first]
    return best


def _stable_surrogate_seed(seed: int, surrogate_id: int, attempt: int) -> int:
    payload = f"{seed}:{surrogate_id}:{attempt}".encode("ascii")
    return int.from_bytes(sha256(payload).digest()[:8], "little", signed=False) % (2**63 - 1)


def _field_hash(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def generate_continuous_h2_surrogates(
    field: np.ndarray,
    count: int,
    seed: int,
    iterations: int,
    h2_offsets: Sequence[Sequence[int]],
    max_attempts_per_id: int,
    swap_iterations: int,
    valid_mask: np.ndarray | None = None,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> tuple[list[np.ndarray], list[dict[str, object]]]:
    """Produce stable-ID continuous surrogates with exact values and H2 gates."""

    source = np.asarray(field, dtype=float)
    if not np.isfinite(source).all():
        raise ValueError("continuous IAAFT surrogates require finite input values")
    if count <= 0 or max_attempts_per_id <= 0:
        raise ValueError("count and max_attempts_per_id must be positive")
    accepted_fields: list[np.ndarray] = []
    audit: list[dict[str, object]] = []
    offsets = _as_offset_array(h2_offsets, source.ndim)
    for surrogate_id in range(1, count + 1):
        accepted = False
        for attempt in range(1, max_attempts_per_id + 1):
            candidate_seed = _stable_surrogate_seed(seed, surrogate_id, attempt)
            candidate = iaaft_surrogate(source, candidate_seed, iterations)
            histogram_exact = bool(np.array_equal(np.sort(source.reshape(-1)), np.sort(candidate.reshape(-1))))
            match = h2_match_diagnostics(
                source, candidate, offsets, valid_mask,
                max_absolute_relative_error=max_absolute_relative_error,
                max_weighted_relative_rmse=max_weighted_relative_rmse,
            )
            # IAAFT candidates that already meet the frozen H2 threshold are
            # retained directly.  Swap optimization is only a repair attempt,
            # never a transformation applied after a passing validation.
            if not bool(match["accepted"]):
                candidate = optimize_h2_by_swaps(
                    source, candidate, offsets, candidate_seed + 1, swap_iterations, valid_mask,
                    max_absolute_relative_error=max_absolute_relative_error,
                    max_weighted_relative_rmse=max_weighted_relative_rmse,
                )
                histogram_exact = bool(np.array_equal(np.sort(source.reshape(-1)), np.sort(candidate.reshape(-1))))
                match = h2_match_diagnostics(
                    source, candidate, offsets, valid_mask,
                    max_absolute_relative_error=max_absolute_relative_error,
                    max_weighted_relative_rmse=max_weighted_relative_rmse,
                )
            candidate_mask = np.isfinite(source) if valid_mask is None else np.asarray(valid_mask, dtype=bool) & np.isfinite(source)
            record: dict[str, object] = {
                "surrogate_id": surrogate_id,
                "generator_type": "continuous_iaaft_swap",
                "seed": candidate_seed,
                "attempt": attempt,
                "input_sha256": _field_hash(source),
                "output_sha256": _field_hash(candidate),
                "histogram_exact": histogram_exact,
                "difference_fraction_from_original": _difference_fraction(source, candidate, candidate_mask),
                "accepted": bool(histogram_exact and match["accepted"]),
                "h2_max_absolute_relative_error": match["maximum_absolute_relative_error"],
                "h2_weighted_relative_rmse": match["weighted_relative_rmse"],
                "h2_pair_count_total": match["pair_count_total"],
                "h2_rows": match["rows"],
            }
            audit.append(record)
            if bool(record["accepted"]):
                accepted_fields.append(candidate)
                accepted = True
                break
        if not accepted:
            last = audit[-1]
            raise RuntimeError(
                f"surrogate_id={surrogate_id} did not satisfy the frozen H2 gate "
                f"(max_error={last['h2_max_absolute_relative_error']}, rmse={last['h2_weighted_relative_rmse']})"
            )
    return accepted_fields, audit


def _strict_majority_from_occupancy(values: np.ndarray, valid_mask: np.ndarray) -> np.ndarray:
    output = np.full(values.shape, np.nan, dtype=float)
    output[valid_mask] = (values[valid_mask] > 0.5).astype(float)
    return output


def _canonical_h2_offsets(offsets: np.ndarray) -> list[tuple[int, int]]:
    """Remove the signed duplicate of each nonzero 2D semivariogram lag."""

    unique: set[tuple[int, int]] = set()
    for raw in np.asarray(offsets, dtype=int):
        offset = (int(raw[0]), int(raw[1]))
        if offset == (0, 0):
            continue
        reverse = (-offset[0], -offset[1])
        unique.add(min(offset, reverse))
    return sorted(unique)


def _pair_indices_for_offset(mask: np.ndarray, offset: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Return flat valid endpoint indices for one signed two-dimensional lag."""

    dy, dx = offset
    ny, nx = mask.shape
    low_y = slice(max(0, -dy), ny - max(0, dy))
    high_y = slice(max(0, dy), ny - max(0, -dy))
    low_x = slice(max(0, -dx), nx - max(0, dx))
    high_x = slice(max(0, dx), nx - max(0, -dx))
    pair_mask = mask[low_y, low_x] & mask[high_y, high_x]
    rows, columns = np.nonzero(pair_mask)
    left_y = rows + max(0, -dy)
    left_x = columns + max(0, -dx)
    right_y = rows + max(0, dy)
    right_x = columns + max(0, dx)
    return (
        (left_y * nx + left_x).astype(np.intp, copy=False),
        (right_y * nx + right_x).astype(np.intp, copy=False),
    )


def _joint_h2_score(
    occupancy_sums: np.ndarray,
    occupancy_target_sums: np.ndarray,
    label_sums: np.ndarray,
    label_target_sums: np.ndarray,
    pair_counts: np.ndarray,
    occupancy_variance: float,
    label_variance: float,
) -> float:
    """Score both discrete H2 constraints on the gate's normalized scale."""

    denominator_occupancy = 2.0 * pair_counts * occupancy_variance
    denominator_label = 2.0 * pair_counts * label_variance
    occupancy_error = np.abs(occupancy_sums - occupancy_target_sums) / denominator_occupancy
    label_error = np.abs(label_sums - label_target_sums) / denominator_label
    weighted_occupancy = float(np.sqrt(np.average(np.square(occupancy_error), weights=pair_counts)))
    weighted_label = float(np.sqrt(np.average(np.square(label_error), weights=pair_counts)))
    # The frozen H2 gate constrains both the worst lag and the pair-weighted
    # RMSE.  Optimizing their sum keeps either field from improving by moving
    # a large mismatch to a single signed lag.
    return float(
        weighted_occupancy
        + weighted_label
        + float(np.max(occupancy_error))
        + float(np.max(label_error))
    )


def _joint_h2_gate(
    occupancy_sums: np.ndarray,
    occupancy_target_sums: np.ndarray,
    label_sums: np.ndarray,
    label_target_sums: np.ndarray,
    pair_counts: np.ndarray,
    occupancy_variance: float,
    label_variance: float,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> bool:
    """Apply the frozen H2 acceptance rule to both discrete representations."""

    occupancy_error = np.abs(occupancy_sums - occupancy_target_sums) / (
        2.0 * pair_counts * occupancy_variance
    )
    label_error = np.abs(label_sums - label_target_sums) / (2.0 * pair_counts * label_variance)
    occupancy_rmse = float(np.sqrt(np.average(np.square(occupancy_error), weights=pair_counts)))
    label_rmse = float(np.sqrt(np.average(np.square(label_error), weights=pair_counts)))
    return bool(
        float(np.max(occupancy_error)) <= max_absolute_relative_error
        and occupancy_rmse <= max_weighted_relative_rmse
        and float(np.max(label_error)) <= max_absolute_relative_error
        and label_rmse <= max_weighted_relative_rmse
    )


def _swap_pair_sum_delta(
    values: np.ndarray,
    labels: np.ndarray,
    first: int,
    second: int,
    pair_indexes: Sequence[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate all semivariogram sum changes caused by a two-cell swap."""

    occupancy_delta = np.zeros(len(pair_indexes), dtype=float)
    label_delta = np.zeros(len(pair_indexes), dtype=float)
    first_value, second_value = float(values[first]), float(values[second])
    first_label, second_label = float(labels[first]), float(labels[second])

    def after(index: int, *, label: bool) -> float:
        if index == first:
            return second_label if label else second_value
        if index == second:
            return first_label if label else first_value
        return float(labels[index]) if label else float(values[index])

    for offset_index, (left, right, left_lookup, right_lookup) in enumerate(pair_indexes):
        # At most four pairs are incident to two cells at one offset.  A set
        # avoids double counting when the two selected cells form a pair.
        incident = {
            int(pair_index)
            for pair_index in (left_lookup[first], right_lookup[first], left_lookup[second], right_lookup[second])
            if pair_index >= 0
        }
        for pair_index in incident:
            left_index, right_index = int(left[pair_index]), int(right[pair_index])
            before_occupancy = (values[right_index] - values[left_index]) ** 2
            after_occupancy = (after(right_index, label=False) - after(left_index, label=False)) ** 2
            occupancy_delta[offset_index] += after_occupancy - before_occupancy
            before_label = (labels[right_index] - labels[left_index]) ** 2
            after_label = (after(right_index, label=True) - after(left_index, label=True)) ** 2
            label_delta[offset_index] += after_label - before_label
    return occupancy_delta, label_delta


def optimize_discrete_h2_by_swaps(
    reference: np.ndarray,
    candidate: np.ndarray,
    h2_offsets: Sequence[Sequence[int]],
    seed: int,
    iterations: int,
    valid_mask: np.ndarray,
    keep_within_h2_gate: bool = False,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> np.ndarray:
    """Anneal an exact discrete multiset against occupancy and label H2 gates.

    A proposal only recalculates pairs touching the two exchanged cells.  The
    resulting cost is proportional to the frozen H2 offset count rather than
    the full field size, which makes the formal 199-field branch practical.
    """

    source = np.asarray(reference, dtype=float)
    current = np.asarray(candidate, dtype=float).copy()
    mask = np.asarray(valid_mask, dtype=bool)
    if source.ndim != 2 or current.shape != source.shape or mask.shape != source.shape:
        raise ValueError("discrete H2 swap optimization requires matching two-dimensional arrays")
    if iterations <= 0:
        return current
    positions = np.flatnonzero(mask.reshape(-1))
    if positions.size < 2:
        return current
    source_flat = source.reshape(-1)
    current_flat = current.reshape(-1)
    source_labels = (source_flat > 0.5).astype(float)
    current_labels = (current_flat > 0.5).astype(float)
    offsets = _canonical_h2_offsets(_as_offset_array(h2_offsets, 2))
    if not offsets:
        return current
    pair_indexes: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for offset in offsets:
        left, right = _pair_indices_for_offset(mask, offset)
        left_lookup = np.full(source_flat.size, -1, dtype=np.intp)
        right_lookup = np.full(source_flat.size, -1, dtype=np.intp)
        pair_ids = np.arange(left.size, dtype=np.intp)
        left_lookup[left] = pair_ids
        right_lookup[right] = pair_ids
        pair_indexes.append((left, right, left_lookup, right_lookup))
    pair_counts = np.asarray([left.size for left, _, _, _ in pair_indexes], dtype=float)
    if np.any(pair_counts <= 0.0):
        raise RuntimeError("a frozen discrete H2 offset has no valid pairs")
    occupancy_variance = float(np.var(source_flat[positions], ddof=0))
    label_variance = float(np.var(source_labels[positions], ddof=0))
    if occupancy_variance <= 0.0 or label_variance <= 0.0:
        raise RuntimeError("discrete H2 fields require positive occupancy and label variance")

    def sums(values: np.ndarray) -> np.ndarray:
        return np.asarray(
            [float(np.sum(np.square(values[right] - values[left]))) for left, right, _, _ in pair_indexes], dtype=float
        )

    occupancy_target = sums(source_flat)
    label_target = sums(source_labels)
    occupancy_current = sums(current_flat)
    label_current = sums(current_labels)
    score = _joint_h2_score(
        occupancy_current, occupancy_target, label_current, label_target, pair_counts,
        occupancy_variance, label_variance,
    )
    currently_within_gate = _joint_h2_gate(
        occupancy_current, occupancy_target, label_current, label_target, pair_counts,
        occupancy_variance, label_variance,
        max_absolute_relative_error=max_absolute_relative_error,
        max_weighted_relative_rmse=max_weighted_relative_rmse,
    )
    if keep_within_h2_gate and not currently_within_gate:
        raise RuntimeError("the discrete H2 walk requires a candidate that already satisfies the frozen H2 gate")
    rng = np.random.default_rng(seed)
    # A bounded temperature is deterministic and deliberately independent of
    # cumulant or connectivity responses.
    for iteration in range(iterations):
        first, second = (int(value) for value in rng.choice(positions, size=2, replace=False))
        if current_flat[first] == current_flat[second]:
            continue
        occupancy_delta, label_delta = _swap_pair_sum_delta(
            current_flat, current_labels, first, second, pair_indexes
        )
        next_occupancy = occupancy_current + occupancy_delta
        next_label = label_current + label_delta
        next_score = _joint_h2_score(
            next_occupancy, occupancy_target, next_label, label_target, pair_counts,
            occupancy_variance, label_variance,
        )
        next_within_gate = _joint_h2_gate(
            next_occupancy, occupancy_target, next_label, label_target, pair_counts,
            occupancy_variance, label_variance,
            max_absolute_relative_error=max_absolute_relative_error,
            max_weighted_relative_rmse=max_weighted_relative_rmse,
        )
        if keep_within_h2_gate:
            # This is a Markov walk on the joint H2-constrained state space.
            # It starts from a valid field and never uses cumulants,
            # connectivity, or any post hoc statistic to choose a proposal.
            accept = next_within_gate
        else:
            temperature = max(1.0e-8, 0.025 * (1.0 - iteration / max(1, iterations)))
            accept = next_score <= score or rng.random() < np.exp((score - next_score) / temperature)
        if accept:
            current_flat[first], current_flat[second] = current_flat[second], current_flat[first]
            current_labels[first], current_labels[second] = current_labels[second], current_labels[first]
            occupancy_current, label_current, score = next_occupancy, next_label, next_score
            currently_within_gate = next_within_gate
    return current


def _difference_fraction(
    reference: np.ndarray,
    candidate: np.ndarray,
    valid_mask: np.ndarray,
) -> float:
    """Return the valid-cell fraction changed by a surrogate proposal."""

    mask = np.asarray(valid_mask, dtype=bool)
    if int(mask.sum()) == 0:
        return 0.0
    return float(np.mean(np.asarray(reference, dtype=float)[mask] != np.asarray(candidate, dtype=float)[mask]))


def generate_discrete_h2_surrogates(
    occupancy: np.ndarray,
    valid_mask: np.ndarray,
    count: int,
    seed: int,
    h2_offsets: Sequence[Sequence[int]],
    max_attempts_per_id: int,
    swap_iterations: int,
    *,
    max_absolute_relative_error: float = 0.10,
    max_weighted_relative_rmse: float = 0.05,
) -> tuple[list[np.ndarray], list[dict[str, object]]]:
    """Generate discrete occupancy fields while matching occupancy and label H2."""

    source = np.asarray(occupancy, dtype=float)
    mask = np.asarray(valid_mask, dtype=bool)
    if source.shape != mask.shape or source.ndim != 2:
        raise ValueError("discrete H2 surrogates require matching 2D occupancy and mask arrays")
    if max_attempts_per_id <= 0:
        raise ValueError("discrete H2 surrogate generation requires max_attempts_per_id > 0")
    if swap_iterations <= 0:
        raise ValueError("discrete H2 surrogate generation requires swap_iterations > 0")
    positions = np.flatnonzero(mask.reshape(-1))
    if positions.size < 2:
        raise ValueError("discrete H2 surrogate field has too few valid cells")
    source_labels = _strict_majority_from_occupancy(source, mask)
    offsets = _as_offset_array(h2_offsets, 2)
    accepted_fields: list[np.ndarray] = []
    audit: list[dict[str, object]] = []
    source_valid_values = source.reshape(-1)[positions].copy()
    for surrogate_id in range(1, count + 1):
        accepted = False
        for attempt in range(1, max_attempts_per_id + 1):
            candidate_seed = _stable_surrogate_seed(seed, surrogate_id, attempt)
            # Start from the exact source field.  The constrained Markov walk
            # below accepts only swaps that retain both frozen H2 conditions.
            # This avoids using translations as either a null model or an
            # initialization shortcut while preserving the discrete histogram.
            candidate = source.copy()
            candidate = optimize_discrete_h2_by_swaps(
                source, candidate, offsets, candidate_seed + 1, swap_iterations, mask, keep_within_h2_gate=True,
                max_absolute_relative_error=max_absolute_relative_error,
                max_weighted_relative_rmse=max_weighted_relative_rmse,
            )
            candidate_labels = _strict_majority_from_occupancy(candidate, mask)
            occupancy_match = h2_match_diagnostics(
                source, candidate, offsets, mask,
                max_absolute_relative_error=max_absolute_relative_error,
                max_weighted_relative_rmse=max_weighted_relative_rmse,
            )
            label_match = h2_match_diagnostics(
                source_labels, candidate_labels, offsets, mask,
                max_absolute_relative_error=max_absolute_relative_error,
                max_weighted_relative_rmse=max_weighted_relative_rmse,
            )
            counts_exact = bool(np.array_equal(np.sort(source_valid_values), np.sort(candidate.reshape(-1)[positions])))
            difference_fraction = _difference_fraction(source, candidate, mask)
            record: dict[str, object] = {
                "surrogate_id": surrogate_id,
                "generator_type": "discrete_h2_constrained_swap_walk",
                "seed": candidate_seed,
                "attempt": attempt,
                "initialization": "source_field",
                "input_sha256": _field_hash(source),
                "output_sha256": _field_hash(candidate),
                "histogram_exact": counts_exact,
                "difference_fraction_from_original": difference_fraction,
                "occupancy_h2_max_absolute_relative_error": occupancy_match["maximum_absolute_relative_error"],
                "occupancy_h2_weighted_relative_rmse": occupancy_match["weighted_relative_rmse"],
                "label_h2_max_absolute_relative_error": label_match["maximum_absolute_relative_error"],
                "label_h2_weighted_relative_rmse": label_match["weighted_relative_rmse"],
                "occupancy_h2_rows": occupancy_match["rows"],
                "label_h2_rows": label_match["rows"],
                "accepted": bool(
                    counts_exact
                    and difference_fraction > 0.0
                    and occupancy_match["accepted"]
                    and label_match["accepted"]
                ),
            }
            audit.append(record)
            if bool(record["accepted"]):
                accepted_fields.append(candidate)
                accepted = True
                break
        if not accepted:
            last = audit[-1]
            raise RuntimeError(
                f"discrete surrogate_id={surrogate_id} did not satisfy occupancy and label H2 gates "
                f"(occupancy_max_error={last['occupancy_h2_max_absolute_relative_error']}, "
                f"occupancy_rmse={last['occupancy_h2_weighted_relative_rmse']}, "
                f"label_max_error={last['label_h2_max_absolute_relative_error']}, "
                f"label_rmse={last['label_h2_weighted_relative_rmse']})"
            )
    return accepted_fields, audit


def surrogate_diversity_diagnostics(fields: Sequence[np.ndarray]) -> dict[str, float | int]:
    """Return frozen diagnostics preventing duplicate fields from filling a run."""

    hashes = [_field_hash(np.asarray(field, dtype=float)) for field in fields]
    pair_rates: list[float] = []
    for left in range(len(fields)):
        for right in range(left + 1, len(fields)):
            a = np.asarray(fields[left])
            b = np.asarray(fields[right])
            pair_rates.append(float(np.mean(a != b)))
    field_count = len(fields)
    duplicate_count = len(hashes) - len(set(hashes))
    return {
        "field_count": field_count,
        "unique_field_count": len(set(hashes)),
        "duplicate_count": duplicate_count,
        "duplicate_fraction": float(duplicate_count / field_count) if field_count else 0.0,
        "minimum_pair_difference_fraction": float(min(pair_rates)) if pair_rates else 1.0,
        "mean_pair_difference_fraction": float(np.mean(pair_rates)) if pair_rates else 1.0,
    }


def symmetric_leave_one_out_test(responses: np.ndarray, numerical_floor: float | None = None) -> dict[str, object]:
    """Permutation-symmetric leave-one-out multivariate surrogate statistic."""

    matrix = np.asarray(responses, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] < 3 or not np.isfinite(matrix).all():
        raise ValueError("responses must be a finite candidates-by-cells matrix with at least three candidates")
    global_range = np.ptp(matrix, axis=0)
    scale = max(1.0, float(np.max(np.abs(matrix))))
    floor = float(numerical_floor) if numerical_floor is not None else np.finfo(float).eps * scale * 64.0
    constant = global_range <= floor
    scores = np.zeros(matrix.shape[0], dtype=float)
    valid_columns = ~constant
    for candidate_index in range(matrix.shape[0]):
        peers = np.delete(matrix, candidate_index, axis=0)
        means = peers.mean(axis=0)
        deviations = peers.std(axis=0, ddof=0)
        if np.any(deviations[valid_columns] <= floor):
            raise RuntimeError("nonconstant response cell has a leave-one-out scale below the frozen numerical floor")
        residual = (matrix[candidate_index, valid_columns] - means[valid_columns]) / deviations[valid_columns]
        scores[candidate_index] = float(np.dot(residual, residual))
    p_value = float(np.count_nonzero(scores >= scores[0]) / scores.size)
    return {
        "scores": scores,
        "original_score": float(scores[0]),
        "empirical_p_value": p_value,
        "numerical_floor": floor,
        "constant_cell_count": int(np.count_nonzero(constant)),
        "active_cell_count": int(np.count_nonzero(valid_columns)),
    }


def holm_step_down(rows: Sequence[dict[str, object]], p_key: str = "p_raw") -> list[dict[str, object]]:
    """Add deterministic Holm adjusted p-values within one already-frozen family."""

    ordered = sorted((dict(row) for row in rows), key=lambda row: (float(row[p_key]), str(row.get("hypothesis_id", ""))))
    family_size = len(ordered)
    running = 0.0
    for rank, row in enumerate(ordered, start=1):
        adjusted = min(1.0, float(row[p_key]) * (family_size - rank + 1))
        running = max(running, adjusted)
        row["holm_rank"] = rank
        row["p_holm"] = running
        row["reject_holm_0_05"] = bool(running < 0.05)
    by_identifier = {str(row.get("hypothesis_id", index)): row for index, row in enumerate(ordered)}
    return [by_identifier[str(row.get("hypothesis_id", index))] for index, row in enumerate(rows)]
