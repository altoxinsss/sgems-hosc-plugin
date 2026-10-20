"""Direct, non-model evaluation for the Section 4.5--4.6 protocol.

The module intentionally contains no classifier, cross-validation, threshold,
or out-of-fold prediction code.  Every result is a direct comparison between
maps at identical moving-window centres and the independently defined channel
facies label.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


COORDINATE_COLUMNS = ("x", "y", "z")
DEFAULT_BASELINES = ("mean", "variance", "high_fraction", "lambda3", "lambda4")
DEFAULT_HIGH_FEATURES = ("kappa4", "lambda5", "lambda6")


@dataclass(frozen=True)
class DirectCase:
    """One coordinate-indexed map family evaluated under a frozen protocol."""

    section: str
    dataset_id: str
    window: str
    weighting: str
    parameter_set: str
    parameter_role: str
    frame: pd.DataFrame
    bootstrap_group: str = ""

    @property
    def identity(self) -> tuple[str, str, str, str, str, str]:
        return (
            self.section,
            self.dataset_id,
            self.window,
            self.weighting,
            self.parameter_set,
            self.parameter_role,
        )


@dataclass(frozen=True)
class BootstrapResult:
    """Point estimates, joint-block intervals, and resampling diagnostics."""

    summary: pd.DataFrame
    diagnostics: dict[str, int | float]


def _require_columns(frame: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{context} is missing required columns: {missing}")


def _ordered(frame: pd.DataFrame) -> pd.DataFrame:
    _require_columns(frame, COORDINATE_COLUMNS, "coordinate frame")
    if frame.duplicated(list(COORDINATE_COLUMNS)).any():
        raise ValueError("Coordinate frame contains duplicate x/y/z locations")
    return frame.sort_values(list(reversed(COORDINATE_COLUMNS))).reset_index(drop=True)


def _assert_same_coordinates(left: pd.DataFrame, right: pd.DataFrame, context: str) -> None:
    left_coordinates = _ordered(left).loc[:, list(COORDINATE_COLUMNS)].to_numpy(dtype=int)
    right_coordinates = _ordered(right).loc[:, list(COORDINATE_COLUMNS)].to_numpy(dtype=int)
    if left_coordinates.shape != right_coordinates.shape or not np.array_equal(left_coordinates, right_coordinates):
        raise ValueError(f"Coordinate support differs for {context}")


def derive_standardized_cumulants(frame: pd.DataFrame, variance_epsilon: float) -> pd.DataFrame:
    """Derive local marginal statistics from C++ ``kappa1`` through ``kappa6``.

    ``lambda4`` is retained solely as the excess-kurtosis representation.  The
    independent ``excess_kurtosis`` expression is evaluated from the recovered
    fourth central moment so that the equality is explicit in the evidence.
    """

    required = [f"kappa{order}" for order in range(1, 7)]
    _require_columns(frame, [*COORDINATE_COLUMNS, *required], "cumulant frame")
    if variance_epsilon <= 0.0:
        raise ValueError("variance_epsilon must be positive")
    result = _ordered(frame).copy()
    variance = pd.to_numeric(result["kappa2"], errors="raise").to_numpy(dtype=float)
    valid = np.isfinite(variance) & (variance > float(variance_epsilon))
    result["mean"] = pd.to_numeric(result["kappa1"], errors="raise").to_numpy(dtype=float)
    result["variance"] = variance
    for order in range(3, 7):
        values = pd.to_numeric(result[f"kappa{order}"], errors="raise").to_numpy(dtype=float)
        standardized = np.full(len(result), np.nan, dtype=float)
        standardized[valid] = values[valid] / variance[valid] ** (order / 2.0)
        result[f"lambda{order}"] = standardized
    recovered_mu4 = pd.to_numeric(result["kappa4"], errors="raise").to_numpy(dtype=float) + 3.0 * variance**2
    excess = np.full(len(result), np.nan, dtype=float)
    excess[valid] = recovered_mu4[valid] / variance[valid] ** 2 - 3.0
    result["excess_kurtosis"] = excess
    result["analysis_variance_valid"] = valid
    return result


def lambda4_equivalence_error(frame: pd.DataFrame) -> dict[str, float | int | bool]:
    """Return the explicit ``lambda4 == excess_kurtosis`` numerical audit."""

    _require_columns(frame, ("lambda4", "excess_kurtosis", "analysis_variance_valid"), "standardized frame")
    mask = frame["analysis_variance_valid"].to_numpy(dtype=bool)
    if not mask.any():
        raise ValueError("No positive-variance rows are available for the lambda4 equivalence audit")
    difference = (
        pd.to_numeric(frame.loc[mask, "lambda4"], errors="raise").to_numpy(dtype=float)
        - pd.to_numeric(frame.loc[mask, "excess_kurtosis"], errors="raise").to_numpy(dtype=float)
    )
    return {
        "valid_count": int(mask.sum()),
        "max_abs_error": float(np.max(np.abs(difference))),
        "mean_abs_error": float(np.mean(np.abs(difference))),
        "finite": bool(np.isfinite(difference).all()),
    }


def build_analysis_support(
    frame: pd.DataFrame,
    *,
    required_features: Sequence[str],
    variance_epsilon: float,
    minimum_ess: float | None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Create a direct-statistics validity mask and auditable exclusion counts.

    Callers attach ``label_valid``, ``cpp_valid`` and
    ``effective_sample_size`` before invoking this function.  The resulting
    ``analysis_valid`` is independent of legacy model feature registries.
    """

    required = [*COORDINATE_COLUMNS, "label_valid", "cpp_valid", "variance", "effective_sample_size", *required_features]
    _require_columns(frame, required, "analysis support frame")
    if variance_epsilon <= 0.0:
        raise ValueError("variance_epsilon must be positive")
    result = _ordered(frame).copy()
    label_valid = result["label_valid"].astype(bool).to_numpy()
    cpp_valid = result["cpp_valid"].astype(bool).to_numpy()
    variance = pd.to_numeric(result["variance"], errors="raise").to_numpy(dtype=float)
    finite_features = np.ones(len(result), dtype=bool)
    for feature in required_features:
        finite_features &= np.isfinite(pd.to_numeric(result[feature], errors="raise").to_numpy(dtype=float))
    variance_valid = np.isfinite(variance) & (variance > float(variance_epsilon))
    ess_valid = np.ones(len(result), dtype=bool)
    if minimum_ess is not None:
        ess = pd.to_numeric(result["effective_sample_size"], errors="raise").to_numpy(dtype=float)
        ess_valid = np.isfinite(ess) & (ess >= float(minimum_ess))
    valid = label_valid & cpp_valid & variance_valid & finite_features & ess_valid
    result["analysis_valid"] = valid
    result["analysis_invalid_reason"] = np.select(
        (
            ~label_valid,
            ~cpp_valid,
            ~variance_valid,
            ~finite_features,
            ~ess_valid,
        ),
        ("facies_tie", "cpp_invalid", "low_variance", "nonfinite_feature", "low_ess"),
        default="included",
    )
    diagnostics = {
        "attempted_count": int(len(result)),
        "included_count": int(valid.sum()),
        "excluded_facies_tie": int((~label_valid).sum()),
        "excluded_cpp_invalid": int((~cpp_valid).sum()),
        "excluded_low_variance": int((~variance_valid).sum()),
        "excluded_nonfinite_feature": int((~finite_features).sum()),
        "excluded_low_ess": int((~ess_valid).sum()),
    }
    return result, diagnostics


def apply_common_support(cases: Sequence[DirectCase]) -> tuple[tuple[DirectCase, ...], dict[str, int]]:
    """Intersect direct-statistics support for paired weighting comparisons."""

    if not cases:
        raise ValueError("At least one case is required for common support")
    reference = _ordered(cases[0].frame)
    _require_columns(reference, ("analysis_valid",), "common-support reference")
    common = reference["analysis_valid"].to_numpy(dtype=bool).copy()
    for case in cases[1:]:
        current = _ordered(case.frame)
        _assert_same_coordinates(reference, current, "common support")
        _require_columns(current, ("analysis_valid",), "common-support case")
        common &= current["analysis_valid"].to_numpy(dtype=bool)
    updated: list[DirectCase] = []
    for case in cases:
        frame = _ordered(case.frame).copy()
        frame["analysis_valid"] = common
        frame["analysis_invalid_reason"] = np.where(common, "included", "not_common_support")
        updated.append(
            DirectCase(
                case.section, case.dataset_id, case.window, case.weighting,
                case.parameter_set, case.parameter_role, frame, case.bootstrap_group,
            )
        )
    return tuple(updated), {"candidate_count": int(len(common)), "common_support_count": int(common.sum())}


def summarize_feature_by_facies(values: np.ndarray, channel: np.ndarray) -> dict[str, float | int]:
    """Return direct group summaries for one channel/non-channel feature map."""

    feature = np.asarray(values, dtype=float)
    label = np.asarray(channel, dtype=bool)
    if feature.ndim != 1 or label.ndim != 1 or feature.size != label.size:
        raise ValueError("Feature and channel labels must be aligned one-dimensional arrays")
    if not np.isfinite(feature).all():
        raise ValueError("Feature values must be finite on analysis support")
    positive = feature[label]
    negative = feature[~label]
    if len(positive) < 2 or len(negative) < 2:
        raise ValueError("Both facies groups require at least two observations")
    return {
        "channel_n": int(len(positive)),
        "nonchannel_n": int(len(negative)),
        "channel_mean": float(np.mean(positive)),
        "nonchannel_mean": float(np.mean(negative)),
        "channel_median": float(np.median(positive)),
        "nonchannel_median": float(np.median(negative)),
        "channel_q1": float(np.quantile(positive, 0.25)),
        "channel_q3": float(np.quantile(positive, 0.75)),
        "nonchannel_q1": float(np.quantile(negative, 0.25)),
        "nonchannel_q3": float(np.quantile(negative, 0.75)),
    }


def standardized_group_contrast(values: np.ndarray, channel: np.ndarray) -> float:
    """Compute signed channel-minus-nonchannel standardized mean difference."""

    feature = np.asarray(values, dtype=float)
    label = np.asarray(channel, dtype=bool)
    positive = feature[label]
    negative = feature[~label]
    if len(positive) < 2 or len(negative) < 2:
        return float("nan")
    denominator = np.sqrt((np.var(positive, ddof=1) + np.var(negative, ddof=1)) / 2.0)
    if not np.isfinite(denominator) or denominator <= 0.0:
        return float("nan")
    return float((np.mean(positive) - np.mean(negative)) / denominator)


def cliffs_delta(values: np.ndarray, channel: np.ndarray) -> float:
    """Compute channel-minus-nonchannel Cliff's delta in O(n log n) time."""

    feature = np.asarray(values, dtype=float)
    label = np.asarray(channel, dtype=bool)
    if feature.ndim != 1 or label.ndim != 1 or len(feature) != len(label):
        raise ValueError("Feature and channel labels must be aligned one-dimensional arrays")
    if not np.isfinite(feature).all():
        raise ValueError("Cliff's delta requires finite values")
    positive = feature[label]
    negative = feature[~label]
    if not len(positive) or not len(negative):
        return float("nan")
    ranks = rankdata(feature, method="average")
    positive_rank_sum = float(ranks[label].sum())
    u_statistic = positive_rank_sum - len(positive) * (len(positive) + 1) / 2.0
    return float(2.0 * u_statistic / (len(positive) * len(negative)) - 1.0)


def spearman_with_occupancy(values: np.ndarray, occupancy: np.ndarray) -> tuple[float, bool]:
    """Return Spearman rho and a validity flag that handles constant columns."""

    feature = np.asarray(values, dtype=float)
    target = np.asarray(occupancy, dtype=float)
    if feature.ndim != 1 or target.ndim != 1 or len(feature) != len(target):
        raise ValueError("Feature and occupancy must be aligned one-dimensional arrays")
    if not np.isfinite(feature).all() or not np.isfinite(target).all():
        return float("nan"), False
    if np.ptp(feature) == 0.0 or np.ptp(target) == 0.0:
        return float("nan"), False
    rho = float(spearmanr(feature, target).statistic)
    return rho, bool(np.isfinite(rho))


def _case_metadata(case: DirectCase) -> dict[str, str]:
    return {
        "section": case.section,
        "dataset_id": case.dataset_id,
        "window": case.window,
        "weighting": case.weighting,
        "parameter_set": case.parameter_set,
        "parameter_role": case.parameter_role,
    }


def compute_direct_effects(cases: Sequence[DirectCase], features: Sequence[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute direct feature summaries, separation effects, and occupancy rho."""

    effect_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for case in cases:
        _require_columns(case.frame, ("analysis_valid", "channel", "channel_occupancy", *features), "direct case")
        selected = _ordered(case.frame).loc[lambda data: data["analysis_valid"].astype(bool)].copy()
        if selected.empty:
            raise ValueError(f"Direct case has no analysis support: {case.identity}")
        channel = selected["channel"].to_numpy(dtype=bool)
        occupancy = selected["channel_occupancy"].to_numpy(dtype=float)
        for feature_name in features:
            values = pd.to_numeric(selected[feature_name], errors="raise").to_numpy(dtype=float)
            summary = summarize_feature_by_facies(values, channel)
            rho, rho_valid = spearman_with_occupancy(values, occupancy)
            metadata = _case_metadata(case)
            summary_rows.append({**metadata, "feature": feature_name, **summary})
            effect_rows.append(
                {
                    **metadata,
                    "feature": feature_name,
                    "effect_E": standardized_group_contrast(values, channel),
                    "separation_S": abs(standardized_group_contrast(values, channel)),
                    "cliffs_delta": cliffs_delta(values, channel),
                    "occupancy_spearman_rho": rho,
                    "occupancy_spearman_valid": rho_valid,
                    "analysis_n": int(len(selected)),
                    "channel_n": int(channel.sum()),
                    "nonchannel_n": int((~channel).sum()),
                }
            )
    return pd.DataFrame(effect_rows), pd.DataFrame(summary_rows)


def annotate_signed_direction_consistency(
    effects: pd.DataFrame,
    *,
    expected_dataset_count: int = 3,
) -> pd.DataFrame:
    """Record whether signed ``E`` has one direction across the three subvolumes.

    Macro averages of ``S=|E|`` are useful for magnitude, but cannot establish
    a common direction.  The annotation stays on the per-dataset effect table
    so reports can show both facts without treating an absolute macro average
    as directional evidence.
    """

    required = {
        "section", "dataset_id", "window", "weighting", "parameter_set",
        "parameter_role", "feature", "effect_E",
    }
    _require_columns(effects, required, "signed-direction effect table")
    if expected_dataset_count <= 0:
        raise ValueError("expected_dataset_count must be positive")
    result = effects.copy()
    group_columns = [
        "section", "window", "weighting", "parameter_set", "parameter_role", "feature",
    ]
    annotations: dict[tuple[str, ...], tuple[int, str, bool | None]] = {}
    for identity, group in result.groupby(group_columns, sort=True, dropna=False):
        ordered = group.sort_values("dataset_id", kind="stable")
        values = ordered["effect_E"].to_numpy(dtype=float)
        signs = ["+" if value > 0.0 else "-" if value < 0.0 else "0" for value in values]
        complete = len(ordered) == int(expected_dataset_count)
        consistent: bool | None
        if not complete:
            consistent = None
        else:
            consistent = bool(np.all(values > 0.0) or np.all(values < 0.0))
        annotations[tuple(str(value) for value in identity)] = (int(len(ordered)), "".join(signs), consistent)
    keys = [tuple(str(row[column]) for column in group_columns) for _, row in result.iterrows()]
    result["signed_E_dataset_count"] = [annotations[key][0] for key in keys]
    result["signed_E_direction_pattern"] = [annotations[key][1] for key in keys]
    result["signed_E_direction_consistent_across_datasets"] = [annotations[key][2] for key in keys]
    return result


def compute_single_map_advantages(
    effects: pd.DataFrame,
    *,
    simple_baselines: Sequence[str] = DEFAULT_BASELINES,
    high_features: Sequence[str] = DEFAULT_HIGH_FEATURES,
) -> pd.DataFrame:
    """Compute A4/A5/A6 after selecting the strongest baseline per map family."""

    required = {
        "section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role", "feature", "separation_S"
    }
    _require_columns(effects, required, "effect table")
    rows: list[dict[str, Any]] = []
    groups = effects.groupby(
        ["section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role"], sort=True, dropna=False
    )
    name_map = {"kappa4": "A4", "lambda5": "A5", "lambda6": "A6"}
    for identity, group in groups:
        section, dataset_id, window, weighting, parameter_set, parameter_role = identity
        baseline = group.loc[group["feature"].isin(simple_baselines)]
        if set(simple_baselines) - set(baseline["feature"]):
            raise ValueError(f"Direct advantage baseline set is incomplete for {identity}")
        baseline_strength = float(baseline["separation_S"].max())
        strongest = str(baseline.loc[baseline["separation_S"].idxmax(), "feature"])
        for feature in high_features:
            selected = group.loc[group["feature"] == feature]
            if len(selected) != 1:
                raise ValueError(f"Direct advantage feature {feature} is missing or duplicated for {identity}")
            row = selected.iloc[0]
            rows.append(
                {
                    "section": section,
                    "dataset_id": dataset_id,
                    "window": window,
                    "weighting": weighting,
                    "parameter_set": parameter_set,
                    "parameter_role": parameter_role,
                    "feature": feature,
                    "advantage_name": name_map[feature],
                    "advantage_A": float(row["separation_S"]) - baseline_strength,
                    "feature_separation_S": float(row["separation_S"]),
                    "strongest_baseline_feature": strongest,
                    "strongest_baseline_S": baseline_strength,
                }
            )
    return pd.DataFrame(rows)


def compute_weighting_effects(effects: pd.DataFrame, advantages: pd.DataFrame) -> pd.DataFrame:
    """Compute W and G for each nonuniform primary weighting versus uniform."""

    rows: list[dict[str, Any]] = []
    primary_effects = effects.loc[
        (effects["parameter_role"] == "primary") & (effects["section"].astype(str) == "4.6")
    ].copy()
    primary_advantages = advantages.loc[
        (advantages["parameter_role"] == "primary") & (advantages["section"].astype(str) == "4.6")
    ].copy()
    for weighting in sorted(set(primary_effects["weighting"]) - {"uniform"}):
        weighted = primary_effects.loc[primary_effects["weighting"] == weighting]
        uniform = primary_effects.loc[primary_effects["weighting"] == "uniform"]
        merged = weighted.merge(
            uniform,
            on=["dataset_id", "window", "feature"],
            suffixes=("_weighted", "_uniform"),
            validate="one_to_one",
        )
        for row in merged.itertuples(index=False):
            rows.append(
                {
                    "section": "4.6",
                    "dataset_id": row.dataset_id,
                    "window": row.window,
                    "weighting": weighting,
                    "parameter_set": row.parameter_set_weighted,
                    "parameter_role": "primary",
                    "feature": row.feature,
                    "effect_name": "W",
                    "effect_value": float(row.separation_S_weighted - row.separation_S_uniform),
                }
            )
        weighted_advantage = primary_advantages.loc[primary_advantages["weighting"] == weighting]
        uniform_advantage = primary_advantages.loc[primary_advantages["weighting"] == "uniform"]
        merged_advantage = weighted_advantage.merge(
            uniform_advantage,
            on=["dataset_id", "window", "feature", "advantage_name"],
            suffixes=("_weighted", "_uniform"),
            validate="one_to_one",
        )
        for row in merged_advantage.itertuples(index=False):
            rows.append(
                {
                    "section": "4.6",
                    "dataset_id": row.dataset_id,
                    "window": row.window,
                    "weighting": weighting,
                    "parameter_set": row.parameter_set_weighted,
                    "parameter_role": "primary",
                    "feature": row.feature,
                    "effect_name": "G",
                    "effect_value": float(row.advantage_A_weighted - row.advantage_A_uniform),
                }
            )
    return pd.DataFrame(rows)


def compute_parameter_sensitivity_effects(effects: pd.DataFrame) -> pd.DataFrame:
    """Compare every prespecified side-wing setting with its fixed main setting.

    The caller is responsible for applying common support before computing the
    effect table.  The two returned contrasts retain both signed ``E`` and
    magnitude ``S`` differences so parameter sensitivity cannot be inferred
    from unpaired point estimates alone.
    """

    required = {
        "section", "dataset_id", "window", "weighting", "parameter_set",
        "parameter_role", "feature", "effect_E", "separation_S",
    }
    _require_columns(effects, required, "parameter-sensitivity effect table")
    primary = effects.loc[
        (effects["section"].astype(str) == "4.6") & (effects["parameter_role"] == "primary")
    ].copy()
    sensitivity = effects.loc[
        (effects["section"].astype(str) == "4.6") & (effects["parameter_role"] == "sensitivity")
    ].copy()
    if sensitivity.empty:
        return pd.DataFrame(
            columns=[
                "section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role",
                "main_parameter_set", "feature", "effect_name", "effect_value",
            ]
        )
    rows: list[dict[str, Any]] = []
    key_columns = ["dataset_id", "window", "weighting", "feature"]
    for weighting in sorted(sensitivity["weighting"].unique()):
        wing = sensitivity.loc[sensitivity["weighting"] == weighting].copy()
        main = primary.loc[primary["weighting"] == weighting].copy()
        if main.empty:
            raise ValueError(f"No primary setting is available for sensitivity weighting {weighting}")
        if main.duplicated(key_columns).any():
            raise ValueError(f"Sensitivity weighting {weighting} has nonunique primary effect rows")
        merged = wing.merge(main, on=key_columns, suffixes=("_wing", "_main"), validate="many_to_one")
        for row in merged.itertuples(index=False):
            base = {
                "section": "4.6",
                "dataset_id": row.dataset_id,
                "window": row.window,
                "weighting": weighting,
                "parameter_set": row.parameter_set_wing,
                "parameter_role": "sensitivity",
                "main_parameter_set": row.parameter_set_main,
                "feature": row.feature,
            }
            rows.extend(
                (
                    {
                        **base,
                        "effect_name": "E_sidewing_minus_main",
                        "effect_value": float(row.effect_E_wing - row.effect_E_main),
                    },
                    {
                        **base,
                        "effect_name": "S_sidewing_minus_main",
                        "effect_value": float(row.separation_S_wing - row.separation_S_main),
                    },
                )
            )
    return pd.DataFrame(rows)


def build_spatial_blocks(frame: pd.DataFrame, block_shape_xyz: tuple[int, int, int]) -> np.ndarray:
    """Assign a stable 3-D block identifier to each centre coordinate."""

    ordered = _ordered(frame)
    shape = tuple(int(value) for value in block_shape_xyz)
    if len(shape) != 3 or any(value <= 0 for value in shape):
        raise ValueError("block_shape_xyz must contain three positive dimensions")
    coordinates = ordered.loc[:, list(COORDINATE_COLUMNS)].to_numpy(dtype=int)
    origin = coordinates.min(axis=0)
    bins = (coordinates - origin) // np.asarray(shape, dtype=int)
    return np.asarray([f"{x}:{y}:{z}" for x, y, z in bins], dtype=object)


def _block_sufficient_statistics(values: np.ndarray, channel: np.ndarray, block_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique, inverse = np.unique(block_ids, return_inverse=True)
    stats = np.zeros((len(unique), 2, 3), dtype=float)
    for group_index, mask in enumerate((channel, ~channel)):
        selected = np.flatnonzero(mask)
        np.add.at(stats[:, group_index, 0], inverse[selected], 1.0)
        np.add.at(stats[:, group_index, 1], inverse[selected], values[selected])
        np.add.at(stats[:, group_index, 2], inverse[selected], values[selected] ** 2)
    return unique, stats


def _effect_from_sufficient_statistics(statistics: np.ndarray) -> float:
    count_positive, sum_positive, square_positive = statistics[0]
    count_negative, sum_negative, square_negative = statistics[1]
    if count_positive < 2.0 or count_negative < 2.0:
        return float("nan")
    mean_positive = sum_positive / count_positive
    mean_negative = sum_negative / count_negative
    variance_positive = (square_positive - sum_positive**2 / count_positive) / (count_positive - 1.0)
    variance_negative = (square_negative - sum_negative**2 / count_negative) / (count_negative - 1.0)
    denominator = np.sqrt((variance_positive + variance_negative) / 2.0)
    if not np.isfinite(denominator) or denominator <= 0.0:
        return float("nan")
    return float((mean_positive - mean_negative) / denominator)


def _long_metric_rows(
    effects: pd.DataFrame,
    advantages: pd.DataFrame,
    weighting: pd.DataFrame,
    parameter_sensitivity: pd.DataFrame | None = None,
) -> pd.DataFrame:
    common = ["section", "dataset_id", "window", "weighting", "parameter_set", "parameter_role", "feature"]
    records: list[pd.DataFrame] = []
    direct = effects.loc[:, [*common, "effect_E", "separation_S"]].melt(
        id_vars=common, value_vars=("effect_E", "separation_S"), var_name="metric_name", value_name="value"
    )
    direct["metric_family"] = "direct_effect"
    records.append(direct)
    if not advantages.empty:
        advantage = advantages.loc[:, [*common, "advantage_name", "advantage_A"]].rename(
            columns={"advantage_name": "metric_name", "advantage_A": "value"}
        )
        advantage["metric_family"] = "single_map_advantage"
        records.append(advantage)
    if not weighting.empty:
        weight = weighting.loc[:, [*common, "effect_name", "effect_value"]].rename(
            columns={"effect_name": "metric_name", "effect_value": "value"}
        )
        weight["metric_family"] = "weighting_effect"
        records.append(weight)
    if parameter_sensitivity is not None and not parameter_sensitivity.empty:
        sensitivity = parameter_sensitivity.loc[:, [*common, "effect_name", "effect_value"]].rename(
            columns={"effect_name": "metric_name", "effect_value": "value"}
        )
        sensitivity["metric_family"] = "parameter_sensitivity"
        records.append(sensitivity)
    return pd.concat(records, ignore_index=True, sort=False)


def _macro_average(records: pd.DataFrame) -> pd.DataFrame:
    group_columns = ["metric_family", "metric_name", "section", "window", "weighting", "parameter_set", "parameter_role", "feature"]
    grouped = records.groupby(group_columns, sort=True, dropna=False, as_index=False)["value"].mean()
    grouped.insert(3, "dataset_id", "macro_equal_dataset")
    return grouped.loc[:, records.columns]


def _case_key(case: DirectCase) -> tuple[str, str, str, str, str, str]:
    return case.identity


def _sampling_group(case: DirectCase) -> tuple[str, ...]:
    if case.bootstrap_group:
        return ("explicit", case.bootstrap_group)
    if case.section == "4.6" and case.parameter_role == "primary":
        # Uniform, coordinate, and value weighting share one intersection mask
        # only for the Section 4.6 paired comparison.
        return ("4.6", case.dataset_id, case.window, "primary")
    return (case.section, case.dataset_id, case.window, case.parameter_role, case.weighting, case.parameter_set)


def joint_spatial_block_bootstrap(
    cases: Sequence[DirectCase],
    *,
    features: Sequence[str],
    simple_baselines: Sequence[str] = DEFAULT_BASELINES,
    high_features: Sequence[str] = DEFAULT_HIGH_FEATURES,
    block_shape_xyz: tuple[int, int, int] = (8, 13, 6),
    replicates: int,
    seed: int,
    max_attempt_multiplier: int = 100,
    include_advantages: bool = True,
    include_weighting_effects: bool = True,
    include_parameter_sensitivity: bool = False,
) -> BootstrapResult:
    """Jointly bootstrap E, S, A, W, and G from cached 3-D block statistics.

    Every primary weighting at a common dataset/window receives the same block
    selection.  Resampling is performed separately inside each subvolume and
    macro estimates are equal averages over those subvolumes.
    """

    if replicates <= 0:
        raise ValueError("bootstrap replicates must be positive")
    if not cases:
        raise ValueError("At least one direct case is required for bootstrap")
    point_effects, _ = compute_direct_effects(cases, features)
    point_advantages = (
        compute_single_map_advantages(point_effects, simple_baselines=simple_baselines, high_features=high_features)
        if include_advantages else pd.DataFrame()
    )
    point_weighting = (
        compute_weighting_effects(point_effects, point_advantages)
        if include_weighting_effects and not point_advantages.empty else pd.DataFrame()
    )
    point_sensitivity = compute_parameter_sensitivity_effects(point_effects) if include_parameter_sensitivity else pd.DataFrame()
    point_records = _long_metric_rows(point_effects, point_advantages, point_weighting, point_sensitivity)
    point_records = pd.concat((point_records, _macro_average(point_records)), ignore_index=True)

    case_statistics: dict[tuple[tuple[str, str, str, str, str, str], str], np.ndarray] = {}
    group_blocks: dict[tuple[str, ...], np.ndarray] = {}
    group_coordinates: dict[tuple[str, ...], np.ndarray] = {}
    case_metadata: dict[tuple[str, str, str, str, str, str], dict[str, str]] = {}
    for case in cases:
        selected = _ordered(case.frame).loc[lambda data: data["analysis_valid"].astype(bool)].copy()
        _require_columns(selected, ("channel", *features), "bootstrap case")
        block_ids = build_spatial_blocks(selected, block_shape_xyz)
        group = _sampling_group(case)
        coordinates = selected.loc[:, list(COORDINATE_COLUMNS)].to_numpy(dtype=int)
        if group in group_coordinates and not np.array_equal(group_coordinates[group], coordinates):
            raise ValueError(f"Paired bootstrap cases do not share common support: {group}")
        group_coordinates[group] = coordinates
        unique_blocks = np.unique(block_ids)
        if group in group_blocks and not np.array_equal(group_blocks[group], unique_blocks):
            raise ValueError(f"Paired bootstrap cases do not share block identities: {group}")
        group_blocks[group] = unique_blocks
        channel = selected["channel"].to_numpy(dtype=bool)
        key = _case_key(case)
        case_metadata[key] = _case_metadata(case)
        for feature in features:
            values = pd.to_numeric(selected[feature], errors="raise").to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite direct feature on analysis support: {key}/{feature}")
            _, stats = _block_sufficient_statistics(values, channel, block_ids)
            case_statistics[(key, feature)] = stats

    point_columns = [
        "metric_family", "metric_name", "section", "dataset_id", "window", "weighting", "parameter_set",
        "parameter_role", "feature",
    ]
    point_lookup = {
        tuple(str(getattr(row, column)) for column in point_columns): float(row.value)
        for row in point_records.itertuples(index=False)
    }
    samples: dict[tuple[str, ...], list[float]] = {key: [] for key in point_lookup}
    generator = np.random.default_rng(seed)
    accepted = 0
    rejected = 0
    attempts = 0
    max_attempts = int(replicates) * int(max_attempt_multiplier)
    while accepted < replicates:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError("Joint spatial bootstrap could not obtain enough finite paired replicates")
        selections = {
            group: generator.integers(0, len(blocks), size=len(blocks)) for group, blocks in group_blocks.items()
        }
        effect_rows: list[dict[str, Any]] = []
        for case in cases:
            key = _case_key(case)
            selection = selections[_sampling_group(case)]
            metadata = case_metadata[key]
            for feature in features:
                effect = _effect_from_sufficient_statistics(case_statistics[(key, feature)][selection].sum(axis=0))
                effect_rows.append(
                    {
                        **metadata,
                        "feature": feature,
                        "effect_E": effect,
                        "separation_S": abs(effect) if np.isfinite(effect) else np.nan,
                    }
                )
        effects = pd.DataFrame(effect_rows)
        if not np.isfinite(effects[["effect_E", "separation_S"]].to_numpy(dtype=float)).all():
            rejected += 1
            continue
        advantages = (
            compute_single_map_advantages(effects, simple_baselines=simple_baselines, high_features=high_features)
            if include_advantages else pd.DataFrame()
        )
        weighting = (
            compute_weighting_effects(effects, advantages)
            if include_weighting_effects and not advantages.empty else pd.DataFrame()
        )
        sensitivity = compute_parameter_sensitivity_effects(effects) if include_parameter_sensitivity else pd.DataFrame()
        records = _long_metric_rows(effects, advantages, weighting, sensitivity)
        records = pd.concat((records, _macro_average(records)), ignore_index=True)
        lookup = {
            tuple(str(getattr(row, column)) for column in point_columns): float(row.value)
            for row in records.itertuples(index=False)
        }
        if set(lookup) != set(samples) or not np.isfinite(np.asarray(list(lookup.values()), dtype=float)).all():
            rejected += 1
            continue
        for key, value in lookup.items():
            samples[key].append(value)
        accepted += 1

    rows: list[dict[str, Any]] = []
    for key, values in sorted(samples.items()):
        vector = np.asarray(values, dtype=float)
        metadata = dict(zip(point_columns, key, strict=True))
        rows.append(
            {
                **metadata,
                "point_estimate": point_lookup[key],
                "bootstrap_mean": float(np.mean(vector)),
                "ci95_low": float(np.quantile(vector, 0.025)),
                "ci95_high": float(np.quantile(vector, 0.975)),
                "bootstrap_replicates": int(accepted),
                "bootstrap_seed": int(seed),
                "block_shape_xyz": "x".join(str(value) for value in block_shape_xyz),
            }
        )
    return BootstrapResult(
        pd.DataFrame(rows),
        {
            "requested_replicates": int(replicates),
            "accepted_replicates": int(accepted),
            "rejected_replicates": int(rejected),
            "attempts": int(attempts),
            "seed": int(seed),
            "block_count_groups": int(len(group_blocks)),
        },
    )


def correlation_matrix(cases: Sequence[DirectCase], features: Sequence[str]) -> pd.DataFrame:
    """Return per-case Spearman matrices for direct-statistics reporting."""

    rows: list[dict[str, Any]] = []
    for case in cases:
        selected = _ordered(case.frame).loc[lambda data: data["analysis_valid"].astype(bool)]
        for left in features:
            for right in features:
                left_values = pd.to_numeric(selected[left], errors="raise").to_numpy(dtype=float)
                right_values = pd.to_numeric(selected[right], errors="raise").to_numpy(dtype=float)
                if np.ptp(left_values) == 0.0 or np.ptp(right_values) == 0.0:
                    rho, valid = float("nan"), False
                else:
                    rho, valid = float(spearmanr(left_values, right_values).statistic), True
                rows.append({**_case_metadata(case), "feature_left": left, "feature_right": right, "spearman_rho": rho, "valid": valid})
    return pd.DataFrame(rows)
