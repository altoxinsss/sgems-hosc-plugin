"""Generate paper-facing direct-statistics tables and figures from frozen CSVs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Direct report input does not exist: {path}")
    return path


def _map_slice(values: pd.DataFrame, feature: str, z: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns = ["x", "y", feature]
    if "analysis_valid" in values.columns:
        columns.append("analysis_valid")
    subset = values.loc[values["z"] == int(z), columns].copy()
    # Boolean facies maps also receive the analysis-support mask.  A float
    # view permits invalid centres to remain visible as NaN pixels.
    subset[feature] = pd.to_numeric(subset[feature], errors="raise").astype(float)
    if "analysis_valid" in subset.columns:
        valid = subset["analysis_valid"].astype(bool)
        subset.loc[~valid, feature] = np.nan
    if subset[feature].notna().sum() == 0:
        raise ValueError(f"No direct feature values for z={z} and {feature}")
    matrix = subset.loc[:, ["x", "y", feature]].pivot(index="y", columns="x", values=feature).sort_index().sort_index(axis=1)
    return matrix.to_numpy(dtype=float), matrix.columns.to_numpy(dtype=int), matrix.index.to_numpy(dtype=int)


def _add_contract(rows: list[dict[str, str]], artifact: Path, *, sources: list[Path], selector: str, scale: str) -> None:
    rows.append(
        {
            "artifact": artifact.name,
            "artifact_sha256": sha256(artifact),
            "source_csv": ";".join(path.name for path in sources),
            "source_sha256": ";".join(sha256(path) for path in sources),
            "filter": selector,
            "coordinate_and_color_rule": scale,
        }
    )


def _forest(axis, frame: pd.DataFrame, *, title: str, xlabel: str) -> None:
    if frame.empty:
        raise ValueError(f"No rows are available for forest plot: {title}")
    ordered = frame.sort_values(["weighting", "window", "metric_name", "feature", "dataset_id"], kind="stable").reset_index(drop=True)
    labels = [f"{row.weighting} {row.window} {row.metric_name} {row.feature} {row.dataset_id}" for row in ordered.itertuples(index=False)]
    values = ordered["point_estimate"].to_numpy(dtype=float)
    low = ordered["ci95_low"].to_numpy(dtype=float)
    high = ordered["ci95_high"].to_numpy(dtype=float)
    y = np.arange(len(ordered))
    axis.errorbar(values, y, xerr=np.vstack((values - low, high - values)), fmt="o", color="#1f5a94", capsize=3)
    axis.axvline(0.0, color="black", linewidth=0.8)
    axis.set_yticks(y, labels, fontsize=7)
    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.grid(axis="x", color="#d9d9d9", linewidth=0.6)


def _with_bootstrap_intervals(
    table: pd.DataFrame,
    bootstrap: pd.DataFrame,
    *,
    scope: str,
    metric_family: str,
    point_metric_column: str | None = None,
) -> pd.DataFrame:
    """Join the precomputed paired bootstrap intervals to a point table."""

    identity = [
        "section", "dataset_id", "window", "weighting", "parameter_set",
        "parameter_role", "feature",
    ]
    intervals = bootstrap.loc[
        (bootstrap["bootstrap_scope"] == scope)
        & (bootstrap["metric_family"] == metric_family),
        [*identity, "metric_name", "ci95_low", "ci95_high", "bootstrap_replicates", "block_shape_xyz"],
    ].copy()
    if intervals.empty:
        raise ValueError(f"Bootstrap scope {scope!r} has no {metric_family!r} intervals")
    if point_metric_column is not None:
        if point_metric_column not in table.columns:
            raise ValueError(f"Point table is missing metric column {point_metric_column!r}")
        intervals = intervals.rename(columns={"metric_name": point_metric_column})
        return table.merge(
            intervals,
            on=[*identity, point_metric_column],
            how="left",
            validate="one_to_one",
        )
    low = intervals.pivot(index=identity, columns="metric_name", values="ci95_low").add_suffix("_ci95_low")
    high = intervals.pivot(index=identity, columns="metric_name", values="ci95_high").add_suffix("_ci95_high")
    metadata = intervals.groupby(identity, dropna=False, as_index=True)[["bootstrap_replicates", "block_shape_xyz"]].first()
    joined = low.join(high).join(metadata).reset_index()
    return table.merge(joined, on=identity, how="left", validate="one_to_one")


def generate_direct_report(run_dir: Path) -> dict[str, Any]:
    """Read frozen direct outputs and generate registered tables and figures."""

    root = Path(run_dir).resolve()
    tables = root / "tables"
    figures = root / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    effects_path = _require(root / "direct_effects.csv")
    summary_path = _require(root / "direct_feature_summary.csv")
    advantages_path = _require(root / "single_map_advantages.csv")
    weighting_path = _require(root / "weighting_effects.csv")
    sensitivity_path = _require(root / "parameter_sensitivity_contrasts.csv")
    bootstrap_path = _require(root / "bootstrap_summary.csv")
    support_path = _require(root / "support_summary.csv")
    values_path = _require(root / "direct_analysis_values.csv")
    mechanism_path = _require(root / "weighting_mechanism_diagnostics.csv")
    effects = pd.read_csv(effects_path)
    summary = pd.read_csv(summary_path)
    advantages = pd.read_csv(advantages_path)
    weighting = pd.read_csv(weighting_path)
    sensitivity = pd.read_csv(sensitivity_path)
    bootstrap = pd.read_csv(bootstrap_path)
    support = pd.read_csv(support_path)
    values = pd.read_csv(values_path)
    mechanism = pd.read_csv(mechanism_path)
    contract_rows: list[dict[str, str]] = []

    section45_table = effects.loc[
        (effects["section"] == 4.5)
        & (effects["window"] == "5x5x3")
        & effects["feature"].isin(["mean", "variance", "high_fraction", "lambda3", "lambda4", "kappa4", "lambda5", "lambda6"])
    ].copy()
    section45_table = _with_bootstrap_intervals(
        section45_table,
        bootstrap,
        scope="primary",
        metric_family="direct_effect",
    )
    section45_table.to_csv(tables / "section45_direct_effects_primary_window.csv", index=False)
    _add_contract(
        contract_rows,
        tables / "section45_direct_effects_primary_window.csv",
        sources=[effects_path, bootstrap_path],
        selector="section=4.5; window=5x5x3; all datasets; direct features",
        scale="direct effect values; no map color scale",
    )
    section46_table = weighting.loc[
        (weighting["section"].astype(str) == "4.6")
        & (weighting["parameter_role"] == "primary")
    ].copy()
    section46_table = _with_bootstrap_intervals(
        section46_table,
        bootstrap,
        scope="primary",
        metric_family="weighting_effect",
        point_metric_column="effect_name",
    )
    section46_table.to_csv(tables / "section46_weighting_effects_primary_window.csv", index=False)
    _add_contract(
        contract_rows,
        tables / "section46_weighting_effects_primary_window.csv",
        sources=[weighting_path, bootstrap_path],
        selector="section=4.6; primary settings; common support",
        scale="paired W/G effects; no map color scale",
    )
    section46_support = support.loc[
        (support["row_type"] == "support_mask")
        & (support["support_scope"] == "section46_all_three_primary_common")
        & (support["section"].astype(str) == "4.6")
        & (support["parameter_role"] == "primary")
        & support["weighting"].isin(["uniform", "coordinate_gaussian", "soft_high_porosity"])
    ].copy()
    uniqueness = ["dataset_id", "window", "weighting", "parameter_set", "parameter_role"]
    if section46_support.empty or section46_support.duplicated(uniqueness).any():
        raise ValueError("Section 4.6 ESS/N diagnostics must contain one common-support primary row per map family")
    section46_support.to_csv(tables / "section46_support_and_ess.csv", index=False)
    _add_contract(
        contract_rows,
        tables / "section46_support_and_ess.csv",
        sources=[support_path],
        selector="section=4.6; common all-three primary support; one row per map family",
        scale="ESS/N and support rate",
    )
    mechanism.to_csv(tables / "section46_weighting_mechanism_diagnostics.csv", index=False)
    _add_contract(
        contract_rows,
        tables / "section46_weighting_mechanism_diagnostics.csv",
        sources=[mechanism_path],
        selector="coordinate Gaussian central x/y/z profiles; soft upper-tail allocated-weight fraction on common support",
        scale="discrete weights and allocated soft-weight fractions",
    )
    section46_sensitivity = _with_bootstrap_intervals(
        sensitivity,
        bootstrap,
        scope="parameter_sensitivity",
        metric_family="parameter_sensitivity",
        point_metric_column="effect_name",
    )
    section46_sensitivity.to_csv(tables / "section46_parameter_sensitivity_contrasts.csv", index=False)
    _add_contract(
        contract_rows,
        tables / "section46_parameter_sensitivity_contrasts.csv",
        sources=[sensitivity_path, bootstrap_path],
        selector="section=4.6; z45_56; side-wing minus main; paired common support",
        scale="E/S side-wing-minus-main with joint block-bootstrap 95% CI",
    )

    representative = values.loc[
        (values["section"].astype(str) == "4.5")
        & (values["dataset_id"] == "z45_56")
        & (values["window"] == "5x5x3")
        & (values["parameter_role"] == "primary")
        & (values["weighting"] == "uniform")
        & (values["z"] == 50)
    ].copy()
    if representative.empty:
        raise ValueError("Representative Section 4.5 slice z45_56/z=50 is absent")
    fields = ("channel", "kappa4", "lambda5", "lambda6")
    figure, axes = plt.subplots(1, 4, figsize=(12.0, 3.1), constrained_layout=True)
    for axis, field in zip(axes, fields, strict=True):
        image, x, y = _map_slice(representative, field, 50)
        if field == "channel":
            artist = axis.imshow(image, origin="lower", cmap="Greys", vmin=0.0, vmax=1.0, aspect="auto")
        else:
            finite = image[np.isfinite(image)]
            bound = float(np.quantile(np.abs(finite), 0.99)) if finite.size else 1.0
            bound = max(bound, np.finfo(float).eps)
            artist = axis.imshow(image, origin="lower", cmap="coolwarm", vmin=-bound, vmax=bound, aspect="auto")
        axis.set_title(field)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.colorbar(artist, ax=axis, shrink=0.8)
    section45_map = figures / "section45_representative_maps.png"
    figure.savefig(section45_map, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        section45_map,
        sources=[values_path],
        selector="z45_56; uniform; 5x5x3; z=50; channel/kappa4/lambda5/lambda6",
        scale="each cumulant uses symmetric 99th-percentile limits; channel uses [0,1]",
    )

    advantage_bootstrap = bootstrap.loc[
        (bootstrap["metric_family"] == "single_map_advantage")
        & (bootstrap["dataset_id"] == "macro_equal_dataset")
        & (bootstrap["section"].astype(str) == "4.5")
        & (bootstrap["weighting"] == "uniform")
        & bootstrap["bootstrap_scope"].isin(["primary", "scale"])
        & bootstrap["feature"].isin(["kappa4", "lambda5", "lambda6"])
    ].copy()
    figure, axis = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)
    _forest(axis, advantage_bootstrap, title="Section 4.5 single-map advantages", xlabel="A = S(feature) - max S(baseline)")
    advantage_figure = figures / "section45_advantage_forest.png"
    figure.savefig(advantage_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        advantage_figure,
        sources=[bootstrap_path],
        selector="section=4.5; metric_family=single_map_advantage; uniform; macro equal dataset; primary/scale bootstrap",
        scale="x-axis direct separation difference and joint block-bootstrap 95% CI",
    )

    correlation_features = ["mean", "variance", "high_fraction", "lambda3", "lambda4", "kappa4", "lambda5", "lambda6"]
    correlation = representative.loc[:, correlation_features].corr(method="spearman")
    figure, axis = plt.subplots(figsize=(6.0, 5.2), constrained_layout=True)
    image = axis.imshow(correlation.to_numpy(dtype=float), cmap="coolwarm", vmin=-1.0, vmax=1.0)
    axis.set_xticks(range(len(correlation_features)), correlation_features, rotation=45, ha="right")
    axis.set_yticks(range(len(correlation_features)), correlation_features)
    figure.colorbar(image, ax=axis, label="Spearman rho")
    correlation_figure = figures / "section45_feature_correlation.png"
    figure.savefig(correlation_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        correlation_figure,
        sources=[values_path],
        selector="z45_56; uniform; 5x5x3; z=50; direct feature correlation",
        scale="fixed Spearman rho range [-1,1]",
    )

    porosity = np.linspace(0.0, 0.5, 201)
    soft_weight = 0.2 + 0.8 / (1.0 + np.exp(-(porosity - 0.25) / 0.02))
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 6.2), constrained_layout=True)
    axis_profiles, axis_soft, axis_tail, axis_ess = axes.flat
    profile = mechanism.loc[
        mechanism["diagnostic_type"] == "coordinate_gaussian_central_axis_profile"
    ].copy()
    if set(profile["axis"]) != {"x", "y", "z"}:
        raise ValueError("Coordinate Gaussian diagnostics require x/y/z central profiles")
    for axis_name, color in (("x", "#1f5a94"), ("y", "#4a8f63"), ("z", "#7a4b97")):
        current = profile.loc[profile["axis"] == axis_name].sort_values("axis_index")
        axis_profiles.plot(current["normalized_offset"], current["weight"], "o-", label=axis_name, color=color)
    axis_profiles.set(xlabel="normalized centre-axis offset", ylabel="weight", title="Coordinate Gaussian central profiles")
    axis_profiles.legend(title="axis")
    axis_soft.plot(porosity, soft_weight, color="#ad5a20")
    axis_soft.set(xlabel="porosity", ylabel="weight", title="Soft high-porosity rule")
    tail = mechanism.loc[
        mechanism["diagnostic_type"] == "soft_high_porosity_upper_tail_weight_fraction"
    ].copy()
    if tail.empty or tail["tail_weight_fraction"].isna().any():
        raise ValueError("Soft upper-tail weight allocation diagnostics are incomplete")
    tail = tail.sort_values("dataset_id")
    tail_bars = axis_tail.bar(tail["dataset_id"], tail["tail_weight_fraction"], color="#ad5a20")
    axis_tail.set(ylabel="weight share for porosity >= 0.25", title="Soft upper-tail allocation", ylim=(0.0, 1.0))
    axis_tail.tick_params(axis="x", labelrotation=20, labelsize=8)
    for bar, value in zip(tail_bars, tail["tail_weight_fraction"], strict=True):
        axis_tail.text(bar.get_x() + bar.get_width() / 2.0, value, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    ess_support = section46_support.loc[:, ["weighting", "ess_over_n_mean"]].copy()
    ess = ess_support.groupby("weighting", sort=False)["ess_over_n_mean"].mean().reindex(
        ["uniform", "coordinate_gaussian", "soft_high_porosity"]
    )
    if ess.isna().any():
        raise ValueError("Section 4.6 ESS/N diagnostics are incomplete")
    bars = axis_ess.bar(ess.index, ess.to_numpy(dtype=float), color=["#555555", "#1f5a94", "#ad5a20"])
    axis_ess.set(ylabel="mean ESS/N", title="Effective sample size")
    axis_ess.set_ylim(0.0, max(1.0, float(ess.max()) * 1.12))
    axis_ess.tick_params(axis="x", labelrotation=25, labelsize=8)
    for bar, value in zip(bars, ess.to_numpy(dtype=float), strict=True):
        axis_ess.text(bar.get_x() + bar.get_width() / 2.0, value, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    weight_definition = figures / "section46_weight_mechanisms_and_ess.png"
    figure.savefig(weight_definition, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        weight_definition,
        sources=[mechanism_path, support_path],
        selector="primary coordinate x/y/z profiles; soft upper-tail allocation; Section 4.6 common-support ESS/N",
        scale="coordinate sigma=1.0; soft threshold=0.25, width=0.02, strength=0.8",
    )

    map_values = values.loc[
        (values["section"].astype(str) == "4.6")
        & (values["dataset_id"] == "z45_56")
        & (values["window"] == "5x5x3")
        & (values["parameter_role"] == "primary")
        & (values["z"] == 50)
        & values["weighting"].isin(["uniform", "coordinate_gaussian", "soft_high_porosity"])
    ].copy()
    pivoted = {name: frame for name, frame in map_values.groupby("weighting", sort=True)}
    if set(pivoted) != {"uniform", "coordinate_gaussian", "soft_high_porosity"}:
        raise ValueError("Section 4.6 representative maps require the three primary weightings")
    lambda_stack = np.concatenate([frame["lambda6"].dropna().to_numpy(dtype=float) for frame in pivoted.values()])
    bound = max(float(np.quantile(np.abs(lambda_stack), 0.99)), np.finfo(float).eps)
    figure, axes = plt.subplots(1, 5, figsize=(15.0, 3.1), constrained_layout=True)
    images = []
    for axis, weighting_name in zip(axes[:3], ("uniform", "coordinate_gaussian", "soft_high_porosity"), strict=True):
        image, _, _ = _map_slice(pivoted[weighting_name], "lambda6", 50)
        images.append(axis.imshow(image, origin="lower", cmap="coolwarm", vmin=-bound, vmax=bound, aspect="auto"))
        axis.set_title(weighting_name)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    uniform_image, _, _ = _map_slice(pivoted["uniform"], "lambda6", 50)
    for axis, weighting_name in zip(axes[3:], ("coordinate_gaussian", "soft_high_porosity"), strict=True):
        image, _, _ = _map_slice(pivoted[weighting_name], "lambda6", 50)
        difference = image - uniform_image
        images.append(axis.imshow(difference, origin="lower", cmap="coolwarm", vmin=-bound, vmax=bound, aspect="auto"))
        axis.set_title(f"{weighting_name} - uniform")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    figure.colorbar(images[0], ax=axes, shrink=0.78, label="lambda6 / difference")
    lambda6_figure = figures / "section46_lambda6_maps.png"
    figure.savefig(lambda6_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        lambda6_figure,
        sources=[values_path],
        selector="z45_56; 5x5x3; z=50; uniform/coordinate/soft lambda6 and paired differences",
        scale="shared symmetric 99th-percentile lambda6 scale across all five panels",
    )

    weighting_bootstrap = bootstrap.loc[
        (bootstrap["metric_family"] == "weighting_effect")
        & (bootstrap["dataset_id"] == "macro_equal_dataset")
        & (bootstrap["section"].astype(str) == "4.6")
        & (bootstrap["bootstrap_scope"] == "primary")
        & (bootstrap["metric_name"].isin(["W", "G"]))
        & (bootstrap["feature"].isin(["kappa4", "lambda5", "lambda6"]))
    ].copy()
    figure, axis = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)
    _forest(axis, weighting_bootstrap, title="Section 4.6 weighting effects", xlabel="W or G relative to uniform")
    weighting_figure = figures / "section46_weight_effect_forest.png"
    figure.savefig(weighting_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)
    _add_contract(
        contract_rows,
        weighting_figure,
        sources=[bootstrap_path],
        selector="section=4.6; metric_family=weighting_effect; W/G; macro equal dataset; primary bootstrap",
        scale="x-axis paired W/G effect and joint block-bootstrap 95% CI",
    )

    contract = pd.DataFrame(contract_rows)
    contract_path = root / "report_artifact_contract.csv"
    contract.to_csv(contract_path, index=False)
    return {
        "artifact_contract": contract_path,
        "artifact_count": int(len(contract)),
        "tables": [str(path.relative_to(root)) for path in tables.glob("*.csv")],
        "figures": [str(path.relative_to(root)) for path in figures.glob("*.png")],
        "direct_effect_rows": int(len(effects)),
        "summary_rows": int(len(summary)),
        "advantage_rows": int(len(advantages)),
        "weighting_effect_rows": int(len(weighting)),
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Generate direct-statistics results tables and figures from a frozen run directory.")
    parser.add_argument("run_dir", type=Path)
    arguments = parser.parse_args()
    result = generate_direct_report(arguments.run_dir)
    print(json.dumps({key: str(value) if isinstance(value, Path) else value for key, value in result.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
