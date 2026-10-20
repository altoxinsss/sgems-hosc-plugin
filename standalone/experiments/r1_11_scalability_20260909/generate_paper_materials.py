from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


EXPECTED_CASE_IDS = (
    "template_workload_r4_d3_a1-1-1_10x10",
    "template_workload_r4_d3_a1-1-1_15x15",
    "template_order_r4_d3_a1-1-1",
    "template_workload_r4_d3_a1-1-1_25x25",
    "window_r4_w5",
    "window_r4_w9",
    "window_r4_w13",
    "template_order_r5_d3_a2-1-1",
    "template_order_r6_d3_a2-2-1",
    "template_order_r7_d3_a2-2-2",
    "template_order_r8_d3_a3-2-2",
    "template_direction_r6_d2_a3-2",
    "template_direction_r6_d4_a2-1-1-1",
    "template_occurrence_r6_d3_a3-1-1",
    "template_occurrence_r6_d3_a1-3-1",
    "template_occurrence_r6_d3_a1-1-3",
)

TEMPLATE_WORKLOAD_CASES = (
    "template_workload_r4_d3_a1-1-1_10x10",
    "template_workload_r4_d3_a1-1-1_15x15",
    "template_order_r4_d3_a1-1-1",
    "template_workload_r4_d3_a1-1-1_25x25",
)
WINDOW_CASES = ("window_r4_w5", "window_r4_w9", "window_r4_w13")
ORDER_CASES = (
    "template_order_r4_d3_a1-1-1",
    "template_order_r5_d3_a2-1-1",
    "template_order_r6_d3_a2-2-1",
    "template_order_r7_d3_a2-2-2",
    "template_order_r8_d3_a3-2-2",
)
DIRECTION_CASES = (
    "template_direction_r6_d2_a3-2",
    "template_order_r6_d3_a2-2-1",
    "template_direction_r6_d4_a2-1-1-1",
)
OCCURRENCE_CASES = (
    "template_occurrence_r6_d3_a3-1-1",
    "template_occurrence_r6_d3_a1-3-1",
    "template_order_r6_d3_a2-2-1",
    "template_occurrence_r6_d3_a1-1-3",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create r=4--8 manuscript tables, figures, and traceable text from one validated R1-11 run."
    )
    parser.add_argument("--run", type=Path, required=True, help="Validated run directory containing CSV evidence.")
    parser.add_argument("--out", type=Path, required=True, help="New, empty material directory.")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, object], key: str) -> float:
    return float(str(row[key]))


def as_bool(raw: str | bool) -> bool:
    return str(raw).strip().lower() == "true"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_run(run: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    required = ("gate.json", "case_matrix.csv", "raw_runs.csv", "route_summary.csv", "pair_summary.csv")
    for name in required:
        require((run / name).is_file(), f"Missing required evidence file: {run / name}")

    gate = json.loads((run / "gate.json").read_text(encoding="utf-8"))
    require(gate.get("overall_pass") is True, "The input run has not passed its numerical/process gate.")
    require(gate.get("case_count") == len(EXPECTED_CASE_IDS), "The input run has an unexpected case count.")
    require(gate.get("route_count") == 2 * len(EXPECTED_CASE_IDS), "The input run has an unexpected route count.")

    cases = read_csv(run / "case_matrix.csv")
    case_ids = {row["case_id"] for row in cases}
    require(case_ids == set(EXPECTED_CASE_IDS), "The input run does not contain the required r=4--8 matrix.")

    raw_rows = read_csv(run / "raw_runs.csv")
    expected_routes = {(case_id, route) for case_id in EXPECTED_CASE_IDS for route in ("bell", "avl")}
    recorded = Counter(
        (row["case_id"], row["route"])
        for row in raw_rows
        if row["phase"] == "recorded" and row["status"] == "PASS"
    )
    require(set(recorded) == expected_routes, "One or more expected routes has no recorded PASS result.")
    for route_key in expected_routes:
        require(recorded[route_key] == 5, f"Expected five recorded runs for {route_key}; found {recorded[route_key]}.")

    pair_rows = read_csv(run / "pair_summary.csv")
    require({row["case_id"] for row in pair_rows} == set(EXPECTED_CASE_IDS), "Pair summary is incomplete.")
    for row in pair_rows:
        require(as_bool(row["pass"]), f"Numerical gate failed for {row['case_id']}.")
        require(int(row.get("nonfinite_value_count", "0")) == 0, f"Non-finite value recorded for {row['case_id']}.")

    route_rows = read_csv(run / "route_summary.csv")
    require({(row["case_id"], row["route"]) for row in route_rows} == expected_routes, "Route summary is incomplete.")
    return pair_rows, route_rows, cases


def enrich_rows(pair_rows: list[dict[str, str]], route_rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    route_index = {(row["case_id"], row["route"]): row for row in route_rows}
    result: dict[str, dict[str, object]] = {}
    for pair in pair_rows:
        case_id = pair["case_id"]
        bell = route_index[(case_id, "bell")]
        avl = route_index[(case_id, "avl")]
        result[case_id] = {
            **pair,
            "bell_sd": float(bell["wall_time_s_std"]),
            "avl_sd": float(avl["wall_time_s_std"]),
            "bell_recorded_runs": int(bell["recorded_runs"]),
            "avl_recorded_runs": int(avl["recorded_runs"]),
        }
    return result


def pm(row: dict[str, object], prefix: str) -> str:
    return f"{value(row, f'{prefix}_wall_time_s_mean'):.3f} +/- {value(row, f'{prefix}_sd'):.3f}"


def peak_pair(row: dict[str, object]) -> str:
    return f"{value(row, 'bell_peak_working_set_mb'):.3f} / {value(row, 'avl_peak_working_set_mb'):.3f}"


def paper_table_row(row: dict[str, object], label: str, template_structure: str = "") -> dict[str, object]:
    return {
        "configuration": label,
        "template_structure": template_structure,
        "output_points": int(str(row["shared_count"])),
        "Bell_mean_s": f"{value(row, 'bell_wall_time_s_mean'):.4f}",
        "Bell_sd_s": f"{value(row, 'bell_sd'):.4f}",
        "AVL_mean_s": f"{value(row, 'avl_wall_time_s_mean'):.4f}",
        "AVL_sd_s": f"{value(row, 'avl_sd'):.4f}",
        "speedup_Bell_over_AVL": f"{value(row, 'speedup_bell_over_avl'):.3f}",
        "Bell_peak_working_set_MiB": f"{value(row, 'bell_peak_working_set_mb'):.3f}",
        "AVL_peak_working_set_MiB": f"{value(row, 'avl_peak_working_set_mb'):.3f}",
        "max_abs_error": f"{value(row, 'max_abs_error'):.6e}",
        "relative_L2_error": f"{value(row, 'relative_l2'):.6e}",
    }


TABLE_FIELDS = [
    "configuration",
    "template_structure",
    "output_points",
    "Bell_mean_s",
    "Bell_sd_s",
    "AVL_mean_s",
    "AVL_sd_s",
    "speedup_Bell_over_AVL",
    "Bell_peak_working_set_MiB",
    "AVL_peak_working_set_MiB",
    "max_abs_error",
    "relative_L2_error",
]


def write_markdown_table(path: Path, rows: list[dict[str, object]], include_structure: bool = False) -> None:
    headers = ["Configuration"]
    keys = ["configuration"]
    if include_structure:
        headers.append("Occurrence vector")
        keys.append("template_structure")
    headers += ["Points", "Bell (s)", "AVL (s)", "Speedup", "Bell / AVL peak working set (MiB)"]
    keys += ["output_points", "Bell_mean_s", "AVL_mean_s", "speedup_Bell_over_AVL", "peak"]
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        rendered = []
        for key in keys:
            if key == "Bell_mean_s":
                rendered.append(f"{row[key]} +/- {row['Bell_sd_s']}")
            elif key == "AVL_mean_s":
                rendered.append(f"{row[key]} +/- {row['AVL_sd_s']}")
            elif key == "peak":
                rendered.append(f"{row['Bell_peak_working_set_MiB']} / {row['AVL_peak_working_set_MiB']}")
            else:
                rendered.append(str(row[key]))
        lines.append("| " + " | ".join(rendered) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("matplotlib is required to create paper figures.") from error
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
        }
    )
    return plt


def plot_runtime_lines(ax, labels: list[str], rows: list[dict[str, object]], title: str, log_y: bool = False) -> None:
    positions = list(range(len(rows)))
    bell = [value(row, "bell_wall_time_s_mean") for row in rows]
    avl = [value(row, "avl_wall_time_s_mean") for row in rows]
    bell_sd = [value(row, "bell_sd") for row in rows]
    avl_sd = [value(row, "avl_sd") for row in rows]
    ax.errorbar(positions, bell, yerr=bell_sd, marker="o", capsize=3, linewidth=1.7, color="#D55E00", label="Bell")
    ax.errorbar(positions, avl, yerr=avl_sd, marker="s", capsize=3, linewidth=1.7, color="#0072B2", label="AVL")
    if log_y:
        ax.set_yscale("log")
    ax.set_xticks(positions, labels)
    ax.set_title(title)
    ax.set_ylabel("Mean wall-clock time (s)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)


def plot_speedup(ax, labels: list[str], rows: list[dict[str, object]], title: str) -> None:
    values = [value(row, "speedup_bell_over_avl") for row in rows]
    positions = list(range(len(rows)))
    bars = ax.bar(positions, values, color="#009E73", width=0.62)
    ax.axhline(1.0, color="#666666", linewidth=1.0, linestyle="--")
    ax.set_xticks(positions, labels)
    ax.set_title(title)
    ax.set_ylabel("Bell / AVL runtime ratio")
    ax.grid(axis="y", alpha=0.25)
    for bar, number in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{number:.1f}x", ha="center", va="bottom", fontsize=8)


def plot_grouped_runtime(ax, labels: list[str], rows: list[dict[str, object]], title: str) -> None:
    positions = list(range(len(rows)))
    width = 0.35
    bell = [value(row, "bell_wall_time_s_mean") for row in rows]
    avl = [value(row, "avl_wall_time_s_mean") for row in rows]
    bell_sd = [value(row, "bell_sd") for row in rows]
    avl_sd = [value(row, "avl_sd") for row in rows]
    ax.bar([x - width / 2 for x in positions], bell, width, yerr=bell_sd, capsize=3, color="#D55E00", label="Bell")
    ax.bar([x + width / 2 for x in positions], avl, width, yerr=avl_sd, capsize=3, color="#0072B2", label="AVL")
    ax.set_xticks(positions, labels)
    ax.set_title(title)
    ax.set_ylabel("Mean wall-clock time (s)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)


def write_figures(figures: Path, rows: dict[str, dict[str, object]]) -> None:
    plt = require_matplotlib()

    template_rows = [rows[case_id] for case_id in TEMPLATE_WORKLOAD_CASES]
    window_rows = [rows[case_id] for case_id in WINDOW_CASES]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.1), constrained_layout=True)
    plot_runtime_lines(
        axes[0],
        [str(int(str(row["workload"]))) for row in template_rows],
        template_rows,
        "(a) Fourth-order template output workload",
    )
    axes[0].set_xlabel("Number of output points")
    plot_runtime_lines(
        axes[1],
        [str(row["workload_label"]).split("x")[0] for row in window_rows],
        window_rows,
        "(b) Fourth-order moving-window size",
    )
    axes[1].set_xlabel("Moving-window side length")
    fig.savefig(figures / "Fig_8_fourth_order_workload.png", dpi=300)
    plt.close(fig)

    order_rows = [rows[case_id] for case_id in ORDER_CASES]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.1), constrained_layout=True)
    order_labels = [f"r={row['order_r']}" for row in order_rows]
    plot_runtime_lines(axes[0], order_labels, order_rows, "(a) Fixed-workload order scaling", log_y=True)
    axes[0].set_xlabel("Cumulant order")
    plot_speedup(axes[1], order_labels, order_rows, "(b) Runtime ratio at each order")
    axes[1].set_xlabel("Cumulant order")
    fig.savefig(figures / "Fig_9_order_scaling_r4_r8.png", dpi=300)
    plt.close(fig)

    direction_rows = [rows[case_id] for case_id in DIRECTION_CASES]
    occurrence_rows = [rows[case_id] for case_id in OCCURRENCE_CASES]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4), constrained_layout=True)
    direction_labels = [
        f"d={row['direction_count']}\n({row['occurrence_vector']})" for row in direction_rows
    ]
    occurrence_labels = [f"({row['occurrence_vector']})" for row in occurrence_rows]
    plot_grouped_runtime(axes[0, 0], direction_labels, direction_rows, "(a) Alternative sixth-order template structures")
    plot_speedup(axes[0, 1], direction_labels, direction_rows, "(b) Runtime ratios for those structures")
    plot_grouped_runtime(axes[1, 0], occurrence_labels, occurrence_rows, "(c) Three-direction occurrence allocations")
    plot_speedup(axes[1, 1], occurrence_labels, occurrence_rows, "(d) Runtime ratios for those allocations")
    fig.savefig(figures / "Fig_10_sixth_order_template_structures.png", dpi=300)
    plt.close(fig)


def write_manuscript_materials(path: Path, rows: dict[str, dict[str, object]]) -> None:
    order_rows = [rows[case_id] for case_id in ORDER_CASES]
    direction_rows = [rows[case_id] for case_id in DIRECTION_CASES]
    occurrence_rows = [rows[case_id] for case_id in OCCURRENCE_CASES]
    largest_error = max(value(row, "max_abs_error") for row in rows.values())
    largest_relative_l2 = max(value(row, "relative_l2") for row in rows.values())
    support_count = min(int(str(row["shared_count"])) for row in rows.values())

    order_lines = []
    for row in order_rows:
        order_lines.append(
            f"- Order {row['order_r']}, occurrence vector ({row['occurrence_vector']}): "
            f"Bell {pm(row, 'bell')} s; AVL {pm(row, 'avl')} s; "
            f"Bell/AVL runtime ratio {value(row, 'speedup_bell_over_avl'):.2f}x; "
            f"recorded peak working set {peak_pair(row)} MiB."
        )

    structure_lines = []
    for row in direction_rows:
        structure_lines.append(
            f"- d={row['direction_count']}, occurrence vector ({row['occurrence_vector']}): "
            f"Bell {pm(row, 'bell')} s; AVL {pm(row, 'avl')} s; "
            f"recorded peak working set {peak_pair(row)} MiB."
        )
    occurrence_lines = []
    for row in occurrence_rows:
        occurrence_lines.append(
            f"- ({row['occurrence_vector']}): Bell {pm(row, 'bell')} s; AVL {pm(row, 'avl')} s; "
            f"recorded peak working set {peak_pair(row)} MiB."
        )

    lines = [
        "# Section 4.3 r=4--8 manuscript materials",
        "",
        "## Figure captions",
        "",
        "Figure 8. Fourth-order runtime across template output workloads and moving-window sizes. "
        "Points and error bars show the mean and standard deviation of five recorded runs; one warm-up run was excluded.",
        "",
        "Figure 9. Runtime scaling from fourth- through eighth-order template cumulants at a fixed 20 x 20 output section. "
        "The Bell panel uses a logarithmic vertical axis because the recorded runtimes span more than one order of magnitude.",
        "",
        "Figure 10. Sixth-order template-complexity comparisons. The upper panels compare configurations that differ in both "
        "direction count and occurrence allocation; the lower panels compare three-direction occurrence allocations.",
        "",
        "## Fixed-workload order-scaling paragraph",
        "",
        "After the fourth-order size benchmarks, the 20 x 20 template was used as the common reference configuration for the order-scaling experiment. "
        "The Stanford VI-E z = 62 slice, the three lag directions, the lag counts, spacing, azimuths, fixed output section, and 400 evaluated positions were held constant; only the occurrence vector was changed to obtain fourth- through eighth-order cumulants.",
        "",
        *order_lines,
        "",
        "## Sixth-order template-complexity paragraph",
        "",
        "The two-, three-, and four-direction cases are comparisons among sixth-order template structures rather than isolated direction-count effects, because the occurrence allocation changes together with the number of directions.",
        "",
        *structure_lines,
        "",
        "For three directions, the occurrence-allocation comparison used four sixth-order configurations:",
        "",
        *occurrence_lines,
        "",
        "## Numerical agreement",
        "",
        f"All cases had matching Bell and AVL coordinate support; every compared case had at least {support_count} shared output points. "
        f"The largest absolute discrepancy was {largest_error:.3e}, and the largest relative L2 discrepancy was {largest_relative_l2:.3e}.",
        "",
        "## Figure and table mapping",
        "",
        "- Table 7 and Figure 8: fourth-order template workload and moving-window size results.",
        "- Table 8 and Figure 9: fixed-workload fourth- through eighth-order series.",
        "- Table 9 and Figure 10: alternative sixth-order template structures and recorded peak working sets.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(out: Path) -> None:
    paths = sorted(path for path in out.rglob("*") if path.is_file() and path.name != "sha256_manifest.csv")
    rows = [
        {
            "path": str(path.relative_to(out)).replace("\\", "/"),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in paths
    ]
    write_csv(out / "sha256_manifest.csv", rows, ["path", "size_bytes", "sha256"])


def main() -> None:
    args = parse_args()
    run = args.run.resolve()
    out = args.out.resolve()
    require(run.is_dir(), f"Run directory does not exist: {run}")
    require(not out.exists(), f"Output directory already exists: {out}")

    pair_rows, route_rows, case_rows = validate_run(run)
    rows = enrich_rows(pair_rows, route_rows)
    tables = out / "tables"
    figures = out / "figures"
    out.mkdir(parents=True)
    tables.mkdir()
    figures.mkdir()

    table_7 = [
        paper_table_row(rows[case_id], f"Template {rows[case_id]['workload_label']} output section")
        for case_id in TEMPLATE_WORKLOAD_CASES
    ] + [
        paper_table_row(rows[case_id], f"Moving window {rows[case_id]['workload_label']}")
        for case_id in WINDOW_CASES
    ]
    table_8 = [
        paper_table_row(rows[case_id], f"Order {rows[case_id]['order_r']}", f"({rows[case_id]['occurrence_vector']})")
        for case_id in ORDER_CASES
    ]
    table_9 = [
        paper_table_row(
            rows[case_id],
            f"Directions {rows[case_id]['direction_count']}",
            f"({rows[case_id]['occurrence_vector']})",
        )
        for case_id in DIRECTION_CASES
    ] + [
        paper_table_row(rows[case_id], "Three-direction allocation", f"({rows[case_id]['occurrence_vector']})")
        for case_id in OCCURRENCE_CASES
    ]

    for filename, table_rows, include_structure in (
        ("Table_7_fourth_order_workload", table_7, False),
        ("Table_8_order_scaling_r4_r8", table_8, True),
        ("Table_9_sixth_order_template_structures", table_9, True),
    ):
        write_csv(tables / f"{filename}.csv", table_rows, TABLE_FIELDS)
        write_markdown_table(tables / f"{filename}.md", table_rows, include_structure)

    metadata = {
        "source_run": str(run),
        "source_gate": json.loads((run / "gate.json").read_text(encoding="utf-8")),
        "source_hashes": {
            name: sha256(run / name)
            for name in ("case_matrix.csv", "raw_runs.csv", "route_summary.csv", "pair_summary.csv", "bell_avl_accuracy.csv", "gate.json")
        },
        "expected_case_ids": list(EXPECTED_CASE_IDS),
        "recorded_repeats_per_route": 5,
        "warmups_per_route": 1,
        "case_rows": case_rows,
    }
    (out / "evidence_reference.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_figures(figures, rows)
    write_manuscript_materials(out / "01_section_4_3_r4_r8_manuscript_materials.md", rows)
    write_manifest(out)
    print(out)


if __name__ == "__main__":
    main()
