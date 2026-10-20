from __future__ import annotations

from dataclasses import dataclass, field
from math import comb

import numpy as np

from window_support import WindowBatch


@dataclass(frozen=True)
class LocalStatistics:
    """One consistent local-statistic set for a complete window support."""

    mean: np.ndarray
    variance: np.ndarray
    central_moment4: np.ndarray
    cumulant4: np.ndarray
    excess_kurtosis: np.ndarray
    high_fraction: np.ndarray
    additional: dict[str, np.ndarray] = field(default_factory=dict)

    def as_dict(self) -> dict[str, np.ndarray]:
        result = {
            "mean": self.mean,
            "variance": self.variance,
            "central_moment4": self.central_moment4,
            "cumulant2": self.variance,
            "cumulant4": self.cumulant4,
            "excess_kurtosis": self.excess_kurtosis,
            "high_fraction": self.high_fraction,
        }
        result.update(self.additional)
        return result


def centered_cumulant(centered: np.ndarray, order: int, weights: np.ndarray | None = None) -> np.ndarray:
    """Evaluate an empirical univariate cumulant from centered samples."""
    if order < 1:
        raise ValueError("Cumulant order must be positive")
    x = np.asarray(centered, dtype=float)
    if x.ndim != 2:
        raise ValueError("centered values must have shape (windows, samples)")
    if weights is None:
        moments = [np.ones(x.shape[0]), np.mean(x, axis=1)]
        moments.extend(np.mean(x**power, axis=1) for power in range(2, order + 1))
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != x.shape:
            raise ValueError("Weight and centered-value shapes differ")
        weight_sum = np.sum(w, axis=1)
        safe_sum = np.where(weight_sum > 0.0, weight_sum, 1.0)
        moments = [np.ones(x.shape[0]), np.sum(w * x, axis=1) / safe_sum]
        moments.extend(np.sum(w * x**power, axis=1) / safe_sum for power in range(2, order + 1))
    cumulants: list[np.ndarray] = [np.zeros(x.shape[0]) for _ in range(order + 1)]
    cumulants[1] = moments[1]
    for current_order in range(2, order + 1):
        value = moments[current_order].copy()
        for j in range(1, current_order):
            value -= comb(current_order - 1, j - 1) * cumulants[j] * moments[current_order - j]
        cumulants[current_order] = value
    return cumulants[order]


def local_statistics_from_values(
    values: np.ndarray,
    high_threshold: float,
    variance_epsilon: float,
    weights: np.ndarray | None = None,
    max_order: int = 6,
) -> LocalStatistics:
    """Compute the complete feature set used by the high-order experiments.

    `excess_kurtosis` and `cumulant4 / variance**2` are equivalent
    parameterizations, not independent information sources.
    """
    if not 2 <= max_order <= 6:
        raise ValueError("max_order must lie in [2, 6]")
    x = np.asarray(values, dtype=float)
    if x.ndim != 2:
        raise ValueError("values must have shape (windows, samples)")
    if weights is None:
        normalized_weights = None
        weight_sum = np.full(x.shape[0], float(x.shape[1]))
        safe_sum = weight_sum
        mean = np.mean(x, axis=1)
        high_fraction = np.mean(x >= high_threshold, axis=1)
    else:
        normalized_weights = np.asarray(weights, dtype=float)
        if normalized_weights.shape != x.shape:
            raise ValueError("Weight and value shapes differ")
        weight_sum = np.sum(normalized_weights, axis=1)
        safe_sum = np.where(weight_sum > 0.0, weight_sum, 1.0)
        mean = np.sum(normalized_weights * x, axis=1) / safe_sum
        high_fraction = np.sum(normalized_weights * (x >= high_threshold), axis=1) / safe_sum
    centered = x - mean[:, None]
    if normalized_weights is None:
        variance = np.mean(centered**2, axis=1)
        central_moment4 = np.mean(centered**4, axis=1)
    else:
        variance = np.sum(normalized_weights * centered**2, axis=1) / safe_sum
        central_moment4 = np.sum(normalized_weights * centered**4, axis=1) / safe_sum
    cumulants = {order: centered_cumulant(centered, order, normalized_weights) for order in range(2, max_order + 1)}
    cumulant4 = cumulants[4] if 4 in cumulants else central_moment4 - 3.0 * variance**2
    valid = (weight_sum > 0.0) & np.isfinite(variance) & (variance > variance_epsilon)
    for values_array in (mean, variance, central_moment4, high_fraction, *cumulants.values()):
        values_array[weight_sum <= 0.0] = np.nan
    excess_kurtosis = np.full_like(variance, np.nan)
    excess_kurtosis[valid] = cumulant4[valid] / variance[valid] ** 2
    skewness = np.full_like(variance, np.nan)
    if 3 in cumulants:
        skewness[valid] = cumulants[3][valid] / variance[valid] ** 1.5
    additional: dict[str, np.ndarray] = {
        "cumulant3": cumulants.get(3, np.full_like(variance, np.nan)),
        "skewness": skewness,
        # Exposed for numerical audit only.  The formal model registry uses
        # excess_kurtosis once and never treats this as a second feature.
        "standardized_kappa4": excess_kurtosis.copy(),
        "valid": valid.astype(bool),
    }
    for order in range(5, max_order + 1):
        standardized = np.full_like(variance, np.nan)
        standardized[valid] = cumulants[order][valid] / variance[valid] ** (order / 2.0)
        additional[f"cumulant{order}"] = cumulants[order]
        additional[f"standardized_kappa{order}"] = standardized
    return LocalStatistics(mean, variance, central_moment4, cumulant4, excess_kurtosis, high_fraction, additional)


def compute_local_statistics(
    windows: WindowBatch,
    high_threshold: float,
    variance_epsilon: float,
    cumulant_order: int = 4,
) -> LocalStatistics:
    # Keep cumulant_order for existing callers. The returned feature set is
    # complete so every model compares exactly the same support.
    return local_statistics_from_values(
        windows.values, high_threshold, variance_epsilon, max_order=max(6, cumulant_order)
    )
