"""Classifier-free matching-moment positive controls for the direct protocol."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import comb
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from controlled_distributions import ControlledExperiment
from formal_protocol import local_statistics_from_values


@dataclass(frozen=True)
class DirectControlledConfig:
    random_seeds: tuple[int, ...]
    window_sample_counts: tuple[int, ...]
    windows_per_class: int
    high_threshold: float
    variance_epsilon: float
    bootstrap_replicates: int

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, Any]) -> "DirectControlledConfig":
        seeds = tuple(int(value) for value in mapping["random_seeds"])
        counts = tuple(int(value) for value in mapping["window_sample_counts"])
        windows = int(mapping.get("direct_windows_per_class", int(mapping["train_windows_per_class"]) + int(mapping["test_windows_per_class"])))
        config = cls(
            random_seeds=seeds,
            window_sample_counts=counts,
            windows_per_class=windows,
            high_threshold=float(mapping.get("high_threshold", 0.25)),
            variance_epsilon=float(mapping.get("variance_epsilon", 1.0e-15)),
            bootstrap_replicates=int(mapping.get("bootstrap_replicates", 10_000)),
        )
        if not config.random_seeds or len(set(config.random_seeds)) != len(config.random_seeds):
            raise ValueError("Direct controlled experiments require unique fixed seeds")
        if not config.window_sample_counts or any(value < 2 for value in config.window_sample_counts):
            raise ValueError("Direct controlled window sample counts must be at least two")
        if config.windows_per_class <= 0 or config.bootstrap_replicates <= 0:
            raise ValueError("Direct controlled window and bootstrap counts must be positive")
        return config

    def as_jsonable(self) -> dict[str, Any]:
        value = asdict(self)
        value["random_seeds"] = list(self.random_seeds)
        value["window_sample_counts"] = list(self.window_sample_counts)
        return value


def _rng(seed: int, experiment_index: int, sample_count: int, class_index: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([seed, experiment_index, sample_count, class_index, 20261011]))


def _sample_non_degenerate(distribution, generator: np.random.Generator, count: int, samples: int, epsilon: float) -> tuple[np.ndarray, int]:
    accepted: list[np.ndarray] = []
    remaining = int(count)
    rejected_total = 0
    while remaining:
        draw = distribution.sample_windows(generator, remaining, samples)
        valid = np.var(draw, axis=1) > epsilon
        accepted.append(draw[valid])
        rejected = int((~valid).sum())
        rejected_total += rejected
        remaining = rejected
    return np.vstack(accepted), rejected_total


def _cumulant_columns(windows: np.ndarray, config: DirectControlledConfig) -> dict[str, np.ndarray]:
    statistics = local_statistics_from_values(
        windows,
        high_threshold=config.high_threshold,
        variance_epsilon=config.variance_epsilon,
        max_order=6,
    ).as_dict()
    result = {
        "kappa1": np.asarray(statistics["mean"], dtype=float),
        "kappa2": np.asarray(statistics["variance"], dtype=float),
        "kappa3": np.asarray(statistics["cumulant3"], dtype=float),
        "kappa4": np.asarray(statistics["cumulant4"], dtype=float),
        "kappa5": np.asarray(statistics["cumulant5"], dtype=float),
        "kappa6": np.asarray(statistics["cumulant6"], dtype=float),
        "high_fraction": np.asarray(statistics["high_fraction"], dtype=float),
    }
    return result


def evaluate_direct_controls(
    experiments: Iterable[ControlledExperiment], config: DirectControlledConfig
) -> pd.DataFrame:
    """Draw finite windows and retain estimator bias/RMSE without model fitting."""

    rows: list[dict[str, Any]] = []
    for experiment_index, experiment in enumerate(experiments):
        target_a = experiment.class_a.affine_cumulants(max_order=6)
        target_b = experiment.class_b.affine_cumulants(max_order=6)
        for sample_count in config.window_sample_counts:
            for seed in config.random_seeds:
                class_a, rejected_a = _sample_non_degenerate(
                    experiment.class_a,
                    _rng(seed, experiment_index, sample_count, 0),
                    config.windows_per_class,
                    sample_count,
                    config.variance_epsilon,
                )
                class_b, rejected_b = _sample_non_degenerate(
                    experiment.class_b,
                    _rng(seed, experiment_index, sample_count, 1),
                    config.windows_per_class,
                    sample_count,
                    config.variance_epsilon,
                )
                estimates_a = _cumulant_columns(class_a, config)
                estimates_b = _cumulant_columns(class_b, config)
                for order in range(1, 7):
                    name = f"kappa{order}"
                    values_a = estimates_a[name]
                    values_b = estimates_b[name]
                    rows.append(
                        {
                            "experiment": experiment.name,
                            "differing_order": experiment.differing_order,
                            "matched_through_order": experiment.matched_through_order,
                            "window_samples": int(sample_count),
                            "seed": int(seed),
                            "order": order,
                            "statistic": name,
                            "windows_per_class": config.windows_per_class,
                            "rejected_degenerate_windows": int(rejected_a + rejected_b),
                            "analytic_class_a": float(target_a[order]),
                            "analytic_class_b": float(target_b[order]),
                            "analytic_delta_a_minus_b": float(target_a[order] - target_b[order]),
                            "estimate_class_a_mean": float(np.mean(values_a)),
                            "estimate_class_b_mean": float(np.mean(values_b)),
                            "estimate_delta_a_minus_b": float(np.mean(values_a) - np.mean(values_b)),
                            "bias_class_a": float(np.mean(values_a) - target_a[order]),
                            "bias_class_b": float(np.mean(values_b) - target_b[order]),
                            "rmse_class_a": float(np.sqrt(np.mean((values_a - target_a[order]) ** 2))),
                            "rmse_class_b": float(np.sqrt(np.mean((values_b - target_b[order]) ** 2))),
                            "high_fraction_class_a": float(np.mean(estimates_a["high_fraction"])),
                            "high_fraction_class_b": float(np.mean(estimates_b["high_fraction"])),
                        }
                    )
    return pd.DataFrame(rows)


def _bootstrap_interval(values: np.ndarray, replicates: int, seed: int) -> tuple[float, float]:
    if len(values) == 1:
        return float(values[0]), float(values[0])
    generator = np.random.default_rng(seed)
    draws = generator.integers(0, len(values), size=(replicates, len(values)))
    means = values[draws].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def summarize_direct_controls(results: pd.DataFrame, config: DirectControlledConfig) -> pd.DataFrame:
    """Summarize finite-window estimator differences over the prespecified seeds."""

    required = {"experiment", "window_samples", "order", "estimate_delta_a_minus_b", "analytic_delta_a_minus_b"}
    if not required.issubset(results.columns):
        raise ValueError("Direct controlled results are incomplete")
    rows: list[dict[str, Any]] = []
    for index, (identity, group) in enumerate(results.groupby(["experiment", "window_samples", "order"], sort=True)):
        experiment, samples, order = identity
        deltas = group["estimate_delta_a_minus_b"].to_numpy(dtype=float)
        analytic = float(group["analytic_delta_a_minus_b"].iloc[0])
        low, high = _bootstrap_interval(deltas, config.bootstrap_replicates, 61000 + index)
        target_order = int(group["differing_order"].iloc[0])
        sign = np.sign(analytic)
        target_detected = bool(order == target_order and sign != 0.0 and sign * low > 0.0)
        rows.append(
            {
                "experiment": experiment,
                "window_samples": int(samples),
                "order": int(order),
                "statistic": str(group["statistic"].iloc[0]),
                "seed_count": int(len(group)),
                "analytic_delta_a_minus_b": analytic,
                "estimate_delta_mean": float(np.mean(deltas)),
                "estimate_delta_ci95_low": low,
                "estimate_delta_ci95_high": high,
                "mean_abs_bias_class_a": float(np.mean(np.abs(group["bias_class_a"].to_numpy(dtype=float)))),
                "mean_abs_bias_class_b": float(np.mean(np.abs(group["bias_class_b"].to_numpy(dtype=float)))),
                "mean_rmse_class_a": float(np.mean(group["rmse_class_a"].to_numpy(dtype=float))),
                "mean_rmse_class_b": float(np.mean(group["rmse_class_b"].to_numpy(dtype=float))),
                "matched_population_order": bool(int(order) <= int(group["matched_through_order"].iloc[0])),
                "target_order": target_order,
                "target_order_detected": target_detected,
            }
        )
    return pd.DataFrame(rows)


def direct_control_gates(
    experiments: Iterable[ControlledExperiment], results: pd.DataFrame, summary: pd.DataFrame
) -> dict[str, object]:
    """Gate analytic matching, complete sampling, and target-order detection."""

    experiment_list = tuple(experiments)
    expected_rows = len(experiment_list) * len(results["window_samples"].unique()) * len(results["seed"].unique()) * 6
    analytic = [experiment.analytic_checks() for experiment in experiment_list]
    targets = summary.loc[summary["order"] == summary["target_order"]]
    target_detected = bool(len(targets) == len(experiment_list) * len(results["window_samples"].unique()) and targets["target_order_detected"].all())
    lower = summary.loc[summary["matched_population_order"]].copy()
    return {
        "analytic_moment_checks_pass": bool(all(bool(row["analytic_gate_pass"]) for row in analytic)),
        "expected_result_rows": int(expected_rows),
        "actual_result_rows": int(len(results)),
        "complete_execution_pass": int(expected_rows) == int(len(results)),
        "all_metrics_finite": bool(np.isfinite(results.select_dtypes(include=[np.number]).to_numpy(dtype=float)).all()),
        "target_order_detection_pass": target_detected,
        "matched_lower_order_diagnostic_mean_abs_delta": float(np.mean(np.abs(lower["estimate_delta_mean"]))) if not lower.empty else float("nan"),
        "overall_pass": False,
    }


def finalize_direct_control_gates(gates: dict[str, object]) -> dict[str, object]:
    """Set the direct positive-control status without finite-window model gates."""

    gates["overall_pass"] = bool(
        gates["analytic_moment_checks_pass"]
        and gates["complete_execution_pass"]
        and gates["all_metrics_finite"]
        and gates["target_order_detection_pass"]
    )
    gates["status"] = "VALIDATION_PASS" if bool(gates["overall_pass"]) else "FAIL"
    return gates
