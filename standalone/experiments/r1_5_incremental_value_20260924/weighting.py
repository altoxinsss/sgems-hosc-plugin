from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping

import numpy as np

from statistics import LocalStatistics, local_statistics_from_values
from window_support import WindowBatch


# These are retained only for reproducing historical analyses.  The formal
# R1-5 protocol admits ``FORMAL_WEIGHTING_MODES`` exclusively; callers that
# need a publication run must validate against that set before requesting a
# C++ map.  Keeping the legacy set explicit prevents compatibility branches
# from being mistaken for current experimental schemes.
FORMAL_WEIGHTING_MODES = frozenset({"uniform", "coordinate_gaussian", "soft_high_porosity"})
LEGACY_WEIGHTING_MODES = frozenset({"value_gaussian", "adaptive_upper_tail"})


@dataclass(frozen=True)
class WeightConfig:
    mode: str = "value_gaussian"
    reference_value: float = 0.26
    bandwidth: float = 0.08
    min_weight: float = 0.0
    normalize: bool = True
    enable_coordinate_gaussian: bool | None = None
    enable_fixed_value_gaussian: bool | None = None
    enable_soft_high_porosity: bool | None = None
    enable_adaptive_upper_tail: bool | None = None
    enable_joint_geometric_mean: bool = False
    high_porosity_threshold: float = 0.25
    soft_transition_width: float = 0.02
    soft_weight_strength: float = 1.0
    adaptive_global_threshold: float = 0.25
    adaptive_mean_offset: float = 1.0
    adaptive_transition_scale: float = 0.25
    adaptive_weight_strength: float = 1.0
    adaptive_std_floor: float = 1.0e-6

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> "WeightConfig":
        fields = cls.__dataclass_fields__
        return cls(**{key: value for key, value in values.items() if key in fields})

    def active_mode(self) -> str:
        switches = {
            "coordinate_gaussian": self.enable_coordinate_gaussian,
            "value_gaussian": self.enable_fixed_value_gaussian,
            "soft_high_porosity": self.enable_soft_high_porosity,
            "adaptive_upper_tail": self.enable_adaptive_upper_tail,
        }
        explicit = any(value is not None for value in switches.values())
        if explicit:
            enabled = [name for name, value in switches.items() if value is True]
            if len(enabled) > 1:
                raise ValueError(f"Only one base weighting switch may be enabled; got {enabled}")
            return enabled[0] if enabled else "uniform"
        aliases = {
            "value": "value_gaussian",
            "high_porosity": "soft_high_porosity",
            "adaptive_tail": "adaptive_upper_tail",
            "coordinate": "coordinate_gaussian",
            "spatial_gaussian": "coordinate_gaussian",
        }
        return aliases.get(self.mode, self.mode)

    def validate(self) -> None:
        mode = self.active_mode()
        supported = FORMAL_WEIGHTING_MODES | LEGACY_WEIGHTING_MODES
        if mode not in supported:
            raise ValueError(f"Unsupported weighting mode: {mode}")
        if mode in {"coordinate_gaussian", "value_gaussian"} and self.bandwidth <= 0.0:
            raise ValueError("bandwidth must be positive")
        if mode == "soft_high_porosity" and self.soft_transition_width <= 0.0:
            raise ValueError("soft_transition_width must be positive")
        if not 0.0 <= self.soft_weight_strength <= 1.0:
            raise ValueError("soft_weight_strength must lie in [0, 1]")
        if not 0.0 <= self.adaptive_weight_strength <= 1.0:
            raise ValueError("adaptive_weight_strength must lie in [0, 1]")
        if self.adaptive_transition_scale < 0.0 or self.adaptive_std_floor <= 0.0:
            raise ValueError("adaptive_transition_scale must be non-negative and adaptive_std_floor positive")
        if self.enable_joint_geometric_mean and mode in {"uniform", "coordinate_gaussian"}:
            raise ValueError("joint geometric mean is valid only for value-based schemes")

    def is_legacy(self) -> bool:
        """Return whether this configuration is outside the formal scheme set."""
        return self.active_mode() in LEGACY_WEIGHTING_MODES or bool(self.enable_joint_geometric_mean)

    def scheme_status(self) -> str:
        """Return an explicit provenance label for reports and compatibility code."""
        return "legacy" if self.is_legacy() else "formal"

    def scheme_id(self) -> str:
        public_name = {
            "value_gaussian": "fixed_value_gaussian",
        }.get(self.active_mode(), self.active_mode())
        identifier = f"{public_name}.v1"
        if self.enable_joint_geometric_mean:
            identifier += "+joint_geometric_mean.v1"
        return identifier


@dataclass(frozen=True)
class WeightedStatisticsResult:
    statistics: LocalStatistics
    weight_sum: np.ndarray
    weight_square_sum: np.ndarray
    effective_sample_size: np.ndarray
    positive_weight_count: np.ndarray


def _sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(np.asarray(values, dtype=float), -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(-clipped))


def coordinate_gaussian_window_weights(
    window_size: tuple[int, int, int],
    source_spacing: tuple[float, float, float],
    bandwidth: float,
    min_weight: float = 0.0,
) -> np.ndarray:
    """Return normalized-coordinate Gaussian weights for one 3-D window.

    The formal definition is based on relative window coordinates, not physical
    distance: exp(-0.5 * sum((delta / half_width)**2) / sigma_s**2).
    ``bandwidth`` is sigma_s and is retained as the public configuration name.
    """
    wx, wy, wz = (int(value) for value in window_size)
    del source_spacing
    if min(wx, wy, wz) <= 0 or bandwidth <= 0.0:
        raise ValueError("window dimensions and sigma_s must be positive")
    dx = (np.arange(wx, dtype=float) - 0.5 * (wx - 1)) / max(1.0, 0.5 * (wx - 1))
    dy = (np.arange(wy, dtype=float) - 0.5 * (wy - 1)) / max(1.0, 0.5 * (wy - 1))
    dz = (np.arange(wz, dtype=float) - 0.5 * (wz - 1)) / max(1.0, 0.5 * (wz - 1))
    z_grid, y_grid, x_grid = np.meshgrid(dz, dy, dx, indexing="ij")
    weights = np.exp(-0.5 * (x_grid**2 + y_grid**2 + z_grid**2) / float(bandwidth) ** 2).reshape(-1)
    weights[weights < float(min_weight)] = 0.0
    return weights


def generate_weights(
    values: np.ndarray,
    config: WeightConfig,
    *,
    window_size: tuple[int, int, int] | None = None,
    source_spacing: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> np.ndarray:
    config.validate()
    mode = config.active_mode()
    x = np.asarray(values, dtype=float)
    if mode == "uniform":
        result = np.ones_like(x)
    elif mode == "coordinate_gaussian":
        if x.ndim != 2 or window_size is None:
            raise ValueError("coordinate_gaussian expects a 2-D window matrix and window_size")
        one_window = coordinate_gaussian_window_weights(
            window_size, source_spacing, float(config.bandwidth), float(config.min_weight)
        )
        if x.shape[1] != len(one_window):
            raise ValueError("coordinate_gaussian window_size does not match the sample count")
        return np.broadcast_to(one_window, x.shape).copy()
    elif mode == "value_gaussian":
        bandwidth = max(float(config.bandwidth), np.finfo(float).eps)
        result = np.exp(-0.5 * ((x - float(config.reference_value)) / bandwidth) ** 2)
    elif mode == "soft_high_porosity":
        membership = _sigmoid((x - float(config.high_porosity_threshold)) / float(config.soft_transition_width))
        strength = float(config.soft_weight_strength)
        result = (1.0 - strength) + strength * membership
    elif mode == "adaptive_upper_tail":
        if x.ndim != 2:
            raise ValueError("adaptive_upper_tail expects a 2-D array of windows by samples")
        local_mean = np.mean(x, axis=1, keepdims=True)
        local_std = np.std(x, axis=1, keepdims=True)
        threshold = np.maximum(
            float(config.adaptive_global_threshold),
            local_mean + float(config.adaptive_mean_offset) * local_std,
        )
        scale_std = np.maximum(local_std, float(config.adaptive_std_floor))
        transition = np.maximum(
            float(config.adaptive_transition_scale) * scale_std,
            float(config.adaptive_std_floor),
        )
        membership = _sigmoid((x - threshold) / transition)
        strength = float(config.adaptive_weight_strength)
        result = (1.0 - strength) + strength * membership
    else:
        raise ValueError(f"Unsupported weighting mode: {mode}")
    result[result < float(config.min_weight)] = 0.0
    return result


def stable_control_seed(*parts: object) -> int:
    """Derive a cross-process-stable 64-bit control seed from protocol keys."""
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little", signed=False)


def matched_ess_uniform_weights(
    reference_weights: np.ndarray,
    seed: int = 0,
    *,
    row_keys: list[object] | tuple[object, ...] | None = None,
    repeat: int = 0,
    return_diagnostics: bool = False,
) -> np.ndarray | tuple[np.ndarray, dict[str, np.ndarray]]:
    """Create a repeated, deterministic uniform-subsample ESS control.

    Every row independently randomizes between ``floor(ESS)`` and
    ``ceil(ESS)`` with the fractional ESS as selection probability, so the
    expected retained count equals the target ESS.  ``row_keys`` should carry
    dataset/window/centre identities in formal runs; row position is only a
    development fallback.
    """
    reference = np.asarray(reference_weights, dtype=float)
    if reference.ndim != 2:
        raise ValueError("reference_weights must have shape (windows, samples)")
    sums = reference.sum(axis=1)
    squares = np.square(reference).sum(axis=1)
    ess = np.divide(sums * sums, squares, out=np.zeros_like(sums), where=squares > 0.0)
    if row_keys is not None and len(row_keys) != len(reference):
        raise ValueError("row_keys must contain one stable identity per window")
    result = np.zeros_like(reference)
    actual = np.zeros(len(reference), dtype=float)
    for row, value in enumerate(ess):
        key = row if row_keys is None else row_keys[row]
        generator = np.random.default_rng(stable_control_seed("matched_ess", seed, repeat, key))
        lower = int(np.clip(np.floor(value), 0, reference.shape[1]))
        upper = int(np.clip(np.ceil(value), 0, reference.shape[1]))
        probability_upper = float(np.clip(value - np.floor(value), 0.0, 1.0))
        count = upper if upper > lower and generator.random() < probability_upper else lower
        if count:
            chosen = generator.choice(reference.shape[1], size=count, replace=False)
            result[row, chosen] = 1.0
        actual[row] = float(count)
    if not return_diagnostics:
        return result
    absolute = np.abs(actual - ess)
    return result, {
        "target_ess": ess,
        "actual_ess": actual,
        "absolute_error": absolute,
        "relative_error": np.divide(absolute, ess, out=np.zeros_like(absolute), where=ess > 0.0),
    }


def permuted_same_weight_controls(
    reference_weights: np.ndarray,
    *,
    seed: int,
    repeat: int,
    row_keys: list[object] | tuple[object, ...],
    association_values: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Permute each nonconstant row while preserving its weight multiset.

    ``association_values`` is either normalized distance from the centre or
    porosity.  It is used solely to record the planned pre/post association
    diagnostic; labels and OOF predictions are intentionally absent.
    """
    reference = np.asarray(reference_weights, dtype=float)
    values = np.asarray(association_values, dtype=float)
    if reference.ndim != 2 or values.shape != reference.shape:
        raise ValueError("reference_weights and association_values must have identical 2-D shapes")
    if len(row_keys) != len(reference):
        raise ValueError("row_keys must contain one stable identity per window")
    output = reference.copy()
    # The archived index matrix makes every analysis-control permutation
    # reproducible without expanding it into a multi-million-row CSV.
    permutation_indices = np.broadcast_to(np.arange(reference.shape[1]), reference.shape).copy()
    changed = np.zeros(len(reference), dtype=bool)
    before = np.full(len(reference), np.nan)
    after = np.full(len(reference), np.nan)
    for row, (weights, association, key) in enumerate(zip(reference, values, row_keys)):
        if np.ptp(weights) <= 0.0:
            continue
        generator = np.random.default_rng(stable_control_seed("weight_permutation", seed, repeat, key))
        for _ in range(32):
            permutation = generator.permutation(len(weights))
            candidate = weights[permutation]
            if not np.array_equal(candidate, weights):
                output[row] = candidate
                permutation_indices[row] = permutation
                changed[row] = True
                break
        if not changed[row]:
            raise RuntimeError("Failed to construct a nonidentity permutation for nonconstant weights")
        if np.std(association) > 0.0:
            before[row] = float(np.corrcoef(weights, association)[0, 1])
            after[row] = float(np.corrcoef(output[row], association)[0, 1])
    return output, {
        "permutation_changed": changed,
        "association_correlation_before": before,
        "association_correlation_after": after,
        "weight_sum_before": np.sum(reference, axis=1),
        "weight_sum_after": np.sum(output, axis=1),
        "weight_square_sum_before": np.sum(reference**2, axis=1),
        "weight_square_sum_after": np.sum(output**2, axis=1),
        "permutation_indices": permutation_indices,
    }


def value_gaussian_weights(values: np.ndarray, config: WeightConfig) -> np.ndarray:
    """Backward-compatible alias; dispatches through the selected weighting scheme."""
    return generate_weights(values, config)


def joint_configuration_weights(component_weights: np.ndarray) -> np.ndarray:
    """Geometric mean across configuration components, matching the C++ template option."""
    weights = np.asarray(component_weights, dtype=float)
    if weights.ndim != 2 or weights.shape[1] == 0:
        raise ValueError("component_weights must have shape (configurations, components)")
    safe = np.maximum(weights, np.finfo(float).tiny)
    return np.exp(np.mean(np.log(safe), axis=1))


def compute_weighted_statistics(
    windows: WindowBatch,
    weights: np.ndarray,
    high_threshold: float,
    variance_epsilon: float,
    cumulant_order: int = 4,
) -> WeightedStatisticsResult:
    x = np.asarray(windows.values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.shape != w.shape:
        raise ValueError(f"Weight shape {w.shape} differs from window shape {x.shape}")
    weight_sum = np.sum(w, axis=1)
    positive_count = np.sum(w > 0.0, axis=1)
    sum_squared = np.sum(w**2, axis=1)
    effective_sample_size = np.full_like(weight_sum, np.nan)
    usable = sum_squared > 0.0
    effective_sample_size[usable] = weight_sum[usable] ** 2 / sum_squared[usable]
    stats = local_statistics_from_values(
        x, high_threshold, variance_epsilon, weights=w, max_order=max(6, cumulant_order)
    )
    return WeightedStatisticsResult(stats, weight_sum, sum_squared, effective_sample_size, positive_count)
