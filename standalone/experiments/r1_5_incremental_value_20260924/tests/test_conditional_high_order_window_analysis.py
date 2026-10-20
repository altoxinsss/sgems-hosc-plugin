from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from conditional_high_order_window_analysis import (
    _association_block_aggregates,
    _association_block_bootstrap,
    _association_from_sufficient_statistics,
    _rank_control_frame,
    _rank_projection,
)


def _synthetic_frame() -> tuple[pd.DataFrame, pd.DataFrame]:
    generator = np.random.default_rng(20261012)
    coordinates = [(x, y, z) for z in range(3) for y in range(6) for x in range(8)]
    count = len(coordinates)
    mean = generator.normal(size=count)
    log_variance = 0.4 * mean + generator.normal(scale=0.8, size=count)
    high_fraction = 0.2 * mean + generator.normal(scale=0.5, size=count)
    lambda3 = -0.2 * mean + generator.normal(scale=0.7, size=count)
    lambda4 = 0.3 * log_variance + generator.normal(scale=0.6, size=count)
    lambda5 = 0.3 * lambda3 - 0.2 * lambda4 + generator.normal(scale=0.9, size=count)
    tail_imbalance = 0.65 * lambda5 + 0.2 * lambda3 + generator.normal(scale=0.75, size=count)
    frame = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    frame["dataset_id"] = "synthetic"
    frame["mean"] = mean
    frame["variance"] = np.exp(log_variance)
    frame["high_fraction"] = high_fraction
    frame["lambda3"] = lambda3
    frame["lambda4"] = lambda4
    frame["lambda5"] = lambda5
    frame["lambda6"] = generator.normal(size=count)
    return frame, pd.DataFrame({"tail_imbalance": tail_imbalance})


def test_block_sufficient_statistics_match_direct_rank_residual_association() -> None:
    frame, metrics = _synthetic_frame()
    controls = ("mean", "variance", "high_fraction", "lambda3", "lambda4")
    rank_controls = _rank_control_frame(frame)
    target = _rank_projection(frame, "lambda5", controls, rank_controls=rank_controls)
    aggregate, diagnostics, metric_fits = _association_block_aggregates(
        frame, target, metrics, rank_controls, (3, 3, 2)
    )

    association, slope, _, _ = _association_from_sufficient_statistics(
        aggregate.xtx.sum(axis=0), aggregate.xtz.sum(axis=0), aggregate.ztz.sum(axis=0)
    )
    expected_association = np.corrcoef(target.residual, metric_fits["tail_imbalance"].residual)[0, 1]
    expected_slope = (
        np.dot(target.residual, metric_fits["tail_imbalance"].residual)
        / np.dot(target.residual, target.residual)
    )

    assert association[0] == pytest.approx(expected_association, abs=1.0e-12)
    assert slope[0] == pytest.approx(expected_slope, abs=1.0e-12)
    assert diagnostics.loc[0, "target_max_abs_residual_control_correlation"] < 1.0e-12
    assert diagnostics.loc[0, "metric_max_abs_residual_control_correlation"] < 1.0e-12


def test_block_bootstrap_is_reproducible_and_reports_macro_average() -> None:
    frame, metrics = _synthetic_frame()
    controls = ("mean", "variance", "high_fraction", "lambda3", "lambda4")
    rank_controls = _rank_control_frame(frame)
    target = _rank_projection(frame, "lambda5", controls, rank_controls=rank_controls)
    aggregate, _, _ = _association_block_aggregates(frame, target, metrics, rank_controls, (3, 3, 2))

    first_summary, first_draws, first_diagnostics = _association_block_bootstrap([aggregate], replicates=24, seed=19)
    second_summary, second_draws, second_diagnostics = _association_block_bootstrap([aggregate], replicates=24, seed=19)

    pd.testing.assert_frame_equal(first_summary, second_summary)
    pd.testing.assert_frame_equal(first_draws, second_draws)
    assert first_diagnostics == second_diagnostics == {"accepted": 24, "attempts": 24}
    per_dataset = first_draws.loc[first_draws["dataset_id"] == "synthetic", "conditional_rank_association"].to_numpy()
    macro = first_draws.loc[first_draws["dataset_id"] == "macro_equal_dataset", "conditional_rank_association"].to_numpy()
    assert np.array_equal(per_dataset, macro)


def test_singular_rank_design_is_rejected() -> None:
    with pytest.raises(ValueError, match="singular"):
        _association_from_sufficient_statistics(
            np.asarray([[2.0, 2.0], [2.0, 2.0]]),
            np.ones((2, 2)),
            np.eye(2),
        )
