from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
import platform
import random
import re
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


ATOL = 1.0e-12
RTOL = 1.0e-10
POLL_INTERVAL_S = 0.005
INTERNAL_TIME_RE = re.compile(r"Program takes\s+([0-9.eE+-]+)\s+s")

# These identifiers are deliberately shared between the evidence report and
# the manuscript material generator.  In particular, the r=4 20 x 20 case is
# measured once and serves both the fourth-order workload series and the
# r=4--8 order-extension series.
ORDER_EXTENSION_CASE_IDS = (
    "template_order_r4_d3_a1-1-1",
    "template_order_r5_d3_a2-1-1",
    "template_order_r6_d3_a2-2-1",
    "template_order_r7_d3_a2-2-2",
    "template_order_r8_d3_a3-2-2",
)
FOURTH_ORDER_TEMPLATE_CASE_IDS = (
    "template_workload_r4_d3_a1-1-1_10x10",
    "template_workload_r4_d3_a1-1-1_15x15",
    "template_order_r4_d3_a1-1-1",
    "template_workload_r4_d3_a1-1-1_25x25",
)
FOURTH_ORDER_WINDOW_CASE_IDS = (
    "window_r4_w5",
    "window_r4_w9",
    "window_r4_w13",
)
SIXTH_ORDER_DIRECTION_CASE_IDS = (
    "template_direction_r6_d2_a3-2",
    "template_order_r6_d3_a2-2-1",
    "template_direction_r6_d4_a2-1-1-1",
)
SIXTH_ORDER_OCCURRENCE_CASE_IDS = (
    "template_occurrence_r6_d3_a3-1-1",
    "template_occurrence_r6_d3_a1-3-1",
    "template_order_r6_d3_a2-2-1",
    "template_occurrence_r6_d3_a1-1-3",
)


@dataclass(frozen=True)
class Case:
    case_id: str
    group: str
    mode: str
    order_r: int
    direction_count: int
    occurrence_vector: str
    workload: int
    workload_label: str
    nod: tuple[int, ...] = ()
    nlag: tuple[int, ...] = ()
    lag: tuple[int, ...] = ()
    azm: tuple[int, ...] = ()
    displayed_dirs: tuple[int, ...] = ()
    fixed_positions: tuple[int, ...] = ()
    window_size: int = 0


@dataclass(frozen=True)
class RouteRun:
    case: Case
    route: str
    algorithm: str


class ProcessMemoryCountersEx(ctypes.Structure):
    _fields_ = [
        ("cb", ctypes.c_ulong),
        ("PageFaultCount", ctypes.c_ulong),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
        ("PrivateUsage", ctypes.c_size_t),
    ]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]
    parser = argparse.ArgumentParser(description="Run the R1-11 standalone HOSC scalability matrix.")
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=repo_root / "build" / "Release" if os.name == "nt" else repo_root / "build",
        help="Directory containing template_cli.exe and movewindow_cli.exe.",
    )
    parser.add_argument(
        "--data",
        type=Path,
        required=True,
        help="Stanford VI-E z=62 columnar grid used by the fixed-workload protocol.",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--paper",
        type=Path,
        default=None,
        help="Read-only provenance input. The script only hashes this file.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def join_ints(values: tuple[int, ...]) -> str:
    return ",".join(str(value) for value in values)


def template_case(
    case_id: str,
    group: str,
    nod: tuple[int, ...],
    nlag: tuple[int, ...],
    azm: tuple[int, ...],
    displayed_dirs: tuple[int, ...],
    fixed_positions: tuple[int, ...] = (),
) -> Case:
    direction_count = len(nod)
    order_r = 1 + sum(nod)
    return Case(
        case_id=case_id,
        group=group,
        mode="template",
        order_r=order_r,
        direction_count=direction_count,
        occurrence_vector=join_ints(nod),
        workload=math.prod(nlag[index - 1] for index in displayed_dirs),
        workload_label="x".join(str(nlag[index - 1]) for index in displayed_dirs),
        nod=nod,
        nlag=nlag,
        lag=tuple(3 for _ in nod),
        azm=azm,
        displayed_dirs=displayed_dirs,
        fixed_positions=fixed_positions,
    )


def window_case(order_r: int, size: int) -> Case:
    return Case(
        case_id=f"window_r{order_r}_w{size}",
        group="window_order_size",
        mode="window",
        order_r=order_r,
        direction_count=0,
        occurrence_vector="",
        workload=size * size,
        workload_label=f"{size}x{size}",
        window_size=size,
    )


def build_cases(smoke: bool) -> list[Case]:
    # One coherent batch underlies both existing fourth-order workload tests
    # and the r=4--8 extension. Direction 1 is fixed at lag 1; directions 2
    # and 3 form the square output section. The r=4 20 x 20 case is shared by
    # the workload and order analyses, and the r=6 d3 (2,2,1) case is shared
    # by the order, direction-count, and occurrence-distribution analyses.
    cases = [
        template_case(
            "template_workload_r4_d3_a1-1-1_10x10",
            "template_workload",
            (1, 1, 1),
            (6, 10, 10),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_workload_r4_d3_a1-1-1_15x15",
            "template_workload",
            (1, 1, 1),
            (6, 15, 15),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_order_r4_d3_a1-1-1",
            "template_order_extension",
            (1, 1, 1),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_workload_r4_d3_a1-1-1_25x25",
            "template_workload",
            (1, 1, 1),
            (6, 25, 25),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        window_case(4, 5),
        window_case(4, 9),
        window_case(4, 13),
        template_case(
            "template_order_r5_d3_a2-1-1",
            "template_order_extension",
            (2, 1, 1),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_order_r6_d3_a2-2-1",
            "template_order_extension",
            (2, 2, 1),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_order_r7_d3_a2-2-2",
            "template_order_extension",
            (2, 2, 2),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_order_r8_d3_a3-2-2",
            "template_order_extension",
            (3, 2, 2),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_direction_r6_d2_a3-2",
            "template_direction_count",
            (3, 2),
            (20, 20),
            (180, 90),
            (1, 2),
        ),
        template_case(
            "template_direction_r6_d4_a2-1-1-1",
            "template_direction_count",
            (2, 1, 1, 1),
            (6, 6, 20, 20),
            (180, 90, 0, 45),
            (3, 4),
            (1, 1),
        ),
        template_case(
            "template_occurrence_r6_d3_a3-1-1",
            "template_occurrence_distribution",
            (3, 1, 1),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_occurrence_r6_d3_a1-3-1",
            "template_occurrence_distribution",
            (1, 3, 1),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
        template_case(
            "template_occurrence_r6_d3_a1-1-3",
            "template_occurrence_distribution",
            (1, 1, 3),
            (6, 20, 20),
            (180, 90, 0),
            (2, 3),
            (1,),
        ),
    ]
    if not smoke:
        return cases

    # Smoke mode retains the full 400-point geometry for an r=4/r=8 parity
    # check, but uses one recorded run when invoked with --repeats 1.
    smoke_ids = {
        "template_order_r4_d3_a1-1-1",
        "template_order_r8_d3_a3-2-2",
    }
    return [case for case in cases if case.case_id in smoke_ids]


def route_runs(cases: list[Case]) -> list[RouteRun]:
    result: list[RouteRun] = []
    for case in cases:
        result.append(RouteRun(case, "bell", "Bell algorithm"))
        result.append(RouteRun(case, "avl", "AVL Based algorithm"))
    return result


def template_config(run: RouteRun, data_path: Path, output_dir: Path) -> str:
    case = run.case
    zeros = tuple(0 for _ in case.nod)
    return f"""[general]
calculate_target=Cumulant
calculate_algorithm={run.algorithm}
source_grid={data_path.as_posix()}
source_property=value
coordinate_dimensions=0
target_property={case.case_id}_{run.route}
target_origin=0,0,0
source_origin=0,0,0
source_spacing=1,1,1

[template_shape]
NOD={join_ints(case.nod)}
NLAG={join_ints(case.nlag)}
LAG={join_ints(case.lag)}
LTOL={join_ints(zeros)}
AZM={join_ints(case.azm)}
ATOL={join_ints(zeros)}
BANDWH={join_ints(zeros)}
DIP={join_ints(zeros)}
DTOL={join_ints(zeros)}
BANDWD={join_ints(zeros)}

[target_cross_sections]
displayed_dirs={join_ints(case.displayed_dirs)}
fixed_positions={join_ints(case.fixed_positions)}

[runtime]
output_dir={output_dir.as_posix()}
"""


def window_config(run: RouteRun, data_path: Path, output_dir: Path) -> str:
    case = run.case
    return f"""[general]
calculate_target=Cumulant
calculate_algorithm={run.algorithm}
source_grid={data_path.as_posix()}
source_property=value
coordinate_dimensions=0
target_property={case.case_id}_{run.route}
target_origin=0,0,0
source_origin=0,0,0
source_spacing=1,1,1

[sliding_window]
moment_order={case.order_r}
window_size={case.window_size},{case.window_size},0
move_step=1,1,0
require_full_window=true
require_complete_window=true
require_anchor_present=true
skip_missing_samples=false
output_at_window_center=true
invert_window_output=false
enable_percentile_clipping=false
minimum_samples=1

[weighted_moment]
weight_center_mode=parameter
weight_center_x=0
weight_center_y=0
weight_center_z=0
normalize_weights=true
zero_first_order_moment=false

[runtime]
output_dir={output_dir.as_posix()}
"""


def read_peak_memory(proc: subprocess.Popen[str]) -> tuple[int, int]:
    if os.name != "nt":
        return 0, 0
    counters = ProcessMemoryCountersEx()
    counters.cb = ctypes.sizeof(counters)
    try:
        ok = ctypes.windll.psapi.GetProcessMemoryInfo(
            int(proc._handle), ctypes.byref(counters), counters.cb
        )
    except (AttributeError, OSError, TypeError, ValueError):
        return 0, 0
    if not ok:
        return 0, 0
    return int(counters.PeakWorkingSetSize), int(counters.PrivateUsage)


def execute(
    exe: Path,
    config_path: Path,
    timeout_s: float,
) -> dict:
    started = time.perf_counter()
    proc = subprocess.Popen(
        [str(exe), str(config_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    peak_working_set = 0
    peak_private = 0
    timed_out = False
    while proc.poll() is None:
        working_set, private_bytes = read_peak_memory(proc)
        peak_working_set = max(peak_working_set, working_set)
        peak_private = max(peak_private, private_bytes)
        if time.perf_counter() - started > timeout_s:
            timed_out = True
            proc.kill()
            break
        time.sleep(POLL_INTERVAL_S)
    stdout, stderr = proc.communicate()
    elapsed = time.perf_counter() - started
    working_set, private_bytes = read_peak_memory(proc)
    peak_working_set = max(peak_working_set, working_set)
    peak_private = max(peak_private, private_bytes)
    internal_match = INTERNAL_TIME_RE.search(stdout)
    return {
        "return_code": proc.returncode,
        "timed_out": timed_out,
        "wall_time_s": elapsed,
        "internal_time_s": float(internal_match.group(1)) if internal_match else "",
        "peak_working_set_mb": peak_working_set / (1024.0 * 1024.0),
        "peak_private_mb": peak_private / (1024.0 * 1024.0),
        "stdout": stdout,
        "stderr": stderr,
    }


def output_file(output_dir: Path, target_property: str) -> Path:
    path = output_dir / f"{target_property}.txt"
    if not path.exists():
        raise FileNotFoundError(f"Expected output was not created: {path}")
    return path


def load_values(path: Path) -> dict[tuple[int, int, int], float]:
    values: dict[tuple[int, int, int], float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x = int(float(parts[-4]))
            y = int(float(parts[-3]))
            z = int(float(parts[-2]))
            value = float(parts[-1])
        except ValueError:
            continue
        values[(x, y, z)] = value
    if not values:
        raise ValueError(f"No output values could be parsed from {path}")
    return values


def compare_outputs(bell_path: Path, avl_path: Path) -> dict:
    bell = load_values(bell_path)
    avl = load_values(avl_path)
    bell_coords = set(bell)
    avl_coords = set(avl)
    shared = sorted(bell_coords & avl_coords)
    support_match = bell_coords == avl_coords
    errors = [avl[coord] - bell[coord] for coord in shared]
    abs_errors = [abs(value) for value in errors]
    squared_errors = [value * value for value in errors]
    bell_squared = [bell[coord] * bell[coord] for coord in shared]
    nonfinite_value_count = sum(
        not math.isfinite(bell[coord]) or not math.isfinite(avl[coord])
        for coord in shared
    )
    failures = [
        coord
        for coord in shared
        if (
            not math.isfinite(bell[coord])
            or not math.isfinite(avl[coord])
            or abs(avl[coord] - bell[coord]) > ATOL + RTOL * abs(bell[coord])
        )
    ]
    if nonfinite_value_count:
        max_abs_error = math.inf
        mean_abs_error = math.inf
        rmse = math.inf
        rel_l2 = math.inf
    else:
        max_abs_error = max(abs_errors) if abs_errors else math.inf
        mean_abs_error = statistics.mean(abs_errors) if abs_errors else math.inf
        rmse = math.sqrt(statistics.mean(squared_errors)) if squared_errors else math.inf
        rel_l2 = (
            math.sqrt(sum(squared_errors)) / math.sqrt(sum(bell_squared))
            if sum(bell_squared) > 0.0
            else (0.0 if not any(abs_errors) else math.inf)
        )
    return {
        "bell_count": len(bell),
        "avl_count": len(avl),
        "shared_count": len(shared),
        "support_match": support_match,
        "max_abs_error": max_abs_error,
        "mean_abs_error": mean_abs_error,
        "rmse": rmse,
        "relative_l2": rel_l2,
        "nonfinite_value_count": nonfinite_value_count,
        "tolerance_failures": len(failures),
        "pass": support_match and nonfinite_value_count == 0 and len(failures) == 0 and bool(shared),
    }


def aggregate(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        if row["phase"] != "recorded" or row["status"] != "PASS":
            continue
        grouped.setdefault((row["case_id"], row["route"]), []).append(row)
    summaries: list[dict] = []
    for (case_id, route), group in sorted(grouped.items()):
        wall = [float(row["wall_time_s"]) for row in group]
        memory = [float(row["peak_working_set_mb"]) for row in group]
        private = [float(row["peak_private_mb"]) for row in group]
        internal = [float(row["internal_time_s"]) for row in group if row["internal_time_s"] != ""]
        first = group[0]
        summaries.append(
            {
                "case_id": case_id,
                "group": first["group"],
                "mode": first["mode"],
                "order_r": first["order_r"],
                "direction_count": first["direction_count"],
                "occurrence_vector": first["occurrence_vector"],
                "workload": first["workload"],
                "workload_label": first["workload_label"],
                "route": route,
                "recorded_runs": len(group),
                "wall_time_s_mean": statistics.mean(wall),
                "wall_time_s_std": statistics.stdev(wall) if len(wall) > 1 else 0.0,
                "wall_time_s_median": statistics.median(wall),
                "internal_time_s_mean": statistics.mean(internal) if internal else "",
                "peak_working_set_mb_max": max(memory),
                "peak_private_mb_max": max(private),
            }
        )
    return summaries


def build_pair_summary(summaries: list[dict], comparisons: list[dict]) -> list[dict]:
    by_key = {(row["case_id"], row["route"]): row for row in summaries}
    comparison_by_case = {row["case_id"]: row for row in comparisons}
    pairs: list[dict] = []
    for case_id in sorted(comparison_by_case):
        bell = by_key.get((case_id, "bell"))
        avl = by_key.get((case_id, "avl"))
        comparison = comparison_by_case[case_id]
        if bell is None or avl is None:
            pairs.append({**comparison, "speedup_bell_over_avl": "", "memory_ratio_bell_over_avl": ""})
            continue
        avl_time = float(avl["wall_time_s_mean"])
        avl_memory = float(avl["peak_working_set_mb_max"])
        pairs.append(
            {
                **comparison,
                "group": bell["group"],
                "mode": bell["mode"],
                "order_r": bell["order_r"],
                "direction_count": bell["direction_count"],
                "occurrence_vector": bell["occurrence_vector"],
                "workload": bell["workload"],
                "workload_label": bell["workload_label"],
                "bell_wall_time_s_mean": bell["wall_time_s_mean"],
                "avl_wall_time_s_mean": avl["wall_time_s_mean"],
                "speedup_bell_over_avl": float(bell["wall_time_s_mean"]) / avl_time if avl_time else math.inf,
                "bell_peak_working_set_mb": bell["peak_working_set_mb_max"],
                "avl_peak_working_set_mb": avl["peak_working_set_mb_max"],
                "memory_ratio_bell_over_avl": float(bell["peak_working_set_mb_max"]) / avl_memory if avl_memory else math.inf,
            }
        )
    return pairs


def format_pair_row(pair: dict | None, label: str) -> str:
    if pair is None:
        return f"| {label} | missing | missing | missing | missing | missing |"
    return (
        f"| {label} | {pair['bell_wall_time_s_mean']:.4f} | "
        f"{pair['avl_wall_time_s_mean']:.4f} | {pair['speedup_bell_over_avl']:.4f} | "
        f"{pair['bell_peak_working_set_mb']:.3f} / "
        f"{pair['avl_peak_working_set_mb']:.3f} | {pair['pass']} |"
    )


def write_run_report(run_root: Path, pairs: list[dict], gate: dict) -> None:
    """Write a compact, self-contained summary after all evidence files exist."""
    by_case = {str(pair["case_id"]): pair for pair in pairs}
    valid_pairs = [pair for pair in pairs if bool(pair.get("pass"))]
    max_abs_error = max((float(pair["max_abs_error"]) for pair in valid_pairs), default=math.inf)
    max_relative_l2 = max((float(pair["relative_l2"]) for pair in valid_pairs), default=math.inf)

    lines = [
        "# R1-11 scalability run report",
        "",
        "## Protocol",
        "",
        "- Stanford VI-E z=62 slice; three-direction order series with NLAG=(6,20,20), LAG=(3,3,3), and AZM=(180,90,0).",
        "- Direction 1 is fixed at its first lag; directions 2 and 3 form the 20 x 20 output section (400 points).",
        "- One warm-up and the configured number of recorded runs are executed for each Bell and AVL route.",
        "- Peak memory means the maximum recorded process working set, reported in MiB.",
        "",
        "## Fourth-order workload and moving-window cases",
        "",
        "| Case | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |",
        "|---|---:|---:|---:|---:|:---:|",
    ]
    for case_id in FOURTH_ORDER_TEMPLATE_CASE_IDS:
        pair = by_case.get(case_id)
        label = f"Template r=4 {pair['workload_label']}" if pair else case_id
        lines.append(format_pair_row(pair, label))
    for case_id in FOURTH_ORDER_WINDOW_CASE_IDS:
        pair = by_case.get(case_id)
        label = f"Moving window r=4 {pair['workload_label']}" if pair else case_id
        lines.append(format_pair_row(pair, label))

    lines.extend([
        "",
        "## Fixed-workload order extension",
        "",
        "| Order | Occurrences by direction | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |",
        "|---:|:---:|---:|---:|---:|---:|:---:|",
    ])
    for case_id in ORDER_EXTENSION_CASE_IDS:
        pair = by_case.get(case_id)
        if pair is None:
            lines.append(f"| missing | missing | missing | missing | missing | missing | missing |")
            continue
        lines.append(
            f"| {pair['order_r']} | ({pair['occurrence_vector']}) | "
            f"{pair['bell_wall_time_s_mean']:.4f} | {pair['avl_wall_time_s_mean']:.4f} | "
            f"{pair['speedup_bell_over_avl']:.4f} | {pair['bell_peak_working_set_mb']:.3f} / "
            f"{pair['avl_peak_working_set_mb']:.3f} | {pair['pass']} |"
        )

    lines.extend([
        "",
        "## Sixth-order template structure",
        "",
        "| Study | Directions | Occurrences | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |",
        "|---|---:|:---:|---:|---:|---:|---:|:---:|",
    ])
    for case_id in SIXTH_ORDER_DIRECTION_CASE_IDS:
        pair = by_case.get(case_id)
        if pair is None:
            continue
        lines.append(
            f"| Direction configuration | {pair['direction_count']} | ({pair['occurrence_vector']}) | "
            f"{pair['bell_wall_time_s_mean']:.4f} | {pair['avl_wall_time_s_mean']:.4f} | "
            f"{pair['speedup_bell_over_avl']:.4f} | {pair['bell_peak_working_set_mb']:.3f} / "
            f"{pair['avl_peak_working_set_mb']:.3f} | {pair['pass']} |"
        )
    for case_id in SIXTH_ORDER_OCCURRENCE_CASE_IDS:
        pair = by_case.get(case_id)
        if pair is None:
            continue
        lines.append(
            f"| Occurrence configuration | {pair['direction_count']} | ({pair['occurrence_vector']}) | "
            f"{pair['bell_wall_time_s_mean']:.4f} | {pair['avl_wall_time_s_mean']:.4f} | "
            f"{pair['speedup_bell_over_avl']:.4f} | {pair['bell_peak_working_set_mb']:.3f} / "
            f"{pair['avl_peak_working_set_mb']:.3f} | {pair['pass']} |"
        )

    lines.extend([
        "",
        "## Validity gates",
        "",
        f"- Cases: {gate['case_count']}; routes: {gate['route_count']}; recorded executions: {gate['recorded_run_count']}.",
        f"- Process gate: {gate['all_processes_pass']}; numerical gate: {gate['all_bell_avl_accuracy_checks_pass']}; overall: {gate['overall_pass']}.",
        f"- Largest absolute Bell/AVL discrepancy: {max_abs_error:.17g}.",
        f"- Largest relative L2 discrepancy: {max_relative_l2:.17g}.",
        "",
        "Detailed provenance is recorded in environment.json, case_matrix.csv, raw_runs.csv, ",
        "route_summary.csv, pair_summary.csv, bell_avl_accuracy.csv, gate.json, and sha256_manifest.csv.",
    ])
    (run_root / "RUN_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def source_manifest(repo_root: Path, extra_paths: list[Path]) -> list[dict]:
    paths = [
        repo_root / "CMakeLists.txt",
        repo_root / "SGEMS_COMPATIBILITY.md",
        repo_root / "src" / "template" / "main.cpp",
        repo_root / "src" / "window" / "main.cpp",
        repo_root / "src" / "common" / "grid_io.cpp",
        repo_root / "src" / "common" / "grid_io.h",
        repo_root / "vendor" / "hosc2008_core" / "CumulantCal.cpp",
        repo_root / "vendor" / "hosc2008_core" / "FunctionCall.cpp",
        repo_root / "vendor" / "hosc2008_core" / "Tools.cpp",
        repo_root / "vendor" / "hoscwindow_core" / "CumulantCal.cpp",
        repo_root / "vendor" / "hoscwindow_core" / "FunctionCall.cpp",
        repo_root / "vendor" / "weighted_moment_core" / "WeightedMomentCore.cpp",
        Path(__file__).resolve(),
        *extra_paths,
    ]
    result = []
    for path in paths:
        if path.exists() and path.is_file():
            result.append(
                {
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "modified_local": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                    "sha256": sha256(path),
                }
            )
    return result


def read_cmake_cache(build_dir: Path) -> dict[str, str]:
    """Record the effective build configuration used for paper-facing timings."""
    cache_path = build_dir / "CMakeCache.txt"
    values: dict[str, str] = {}
    if not cache_path.exists():
        return values
    wanted = {
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_FLAGS",
        "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_GENERATOR",
    }
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        lhs, value = line.split("=", 1)
        key = lhs.split(":", 1)[0]
        if key in wanted:
            values[key] = value
    return values


def main() -> int:
    args = parse_args()
    if args.repeats < 1 or args.warmups < 0:
        raise ValueError("repeats must be positive and warmups must be non-negative")

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]
    build_dir = args.build_dir.resolve()
    data_path = args.data.resolve()
    paper_path = args.paper.resolve() if args.paper is not None else None
    suffix = ".exe" if os.name == "nt" else ""
    template_exe = build_dir / ("template_cli" + suffix)
    window_exe = build_dir / ("movewindow_cli" + suffix)
    for required in (template_exe, window_exe, data_path):
        if not required.exists():
            raise FileNotFoundError(required)

    run_stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    run_root = repo_root / "output" / "r1_11_scalability_20260909" / run_stamp
    configs_dir = run_root / "configs"
    logs_dir = run_root / "logs"
    outputs_dir = run_root / "outputs"
    for directory in (configs_dir, logs_dir, outputs_dir):
        directory.mkdir(parents=True, exist_ok=False)

    cases = build_cases(args.smoke)
    routes = route_runs(cases)
    random.Random(args.seed).shuffle(routes)

    environment = {
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": sys.version,
        "repo_root": str(repo_root),
        "build_dir": str(build_dir),
        "cmake_cache": read_cmake_cache(build_dir),
        "data": str(data_path),
        "paper_read_only": str(paper_path) if paper_path is not None else None,
        "repeats": args.repeats,
        "warmups": args.warmups,
        "timeout_s": args.timeout,
        "seed": args.seed,
        "smoke": args.smoke,
        "atol": ATOL,
        "rtol": RTOL,
        "sgems_compatibility": {
            "require_full_window": True,
            "require_complete_window": True,
            "require_anchor_present": True,
            "skip_missing_samples": False,
            "output_at_window_center": True,
            "source_origin": [0, 0, 0],
            "source_spacing": [1, 1, 1],
            "invert_window_output": False,
            "enable_percentile_clipping": False,
        },
    }
    (run_root / "environment.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    case_rows = []
    for case in cases:
        row = asdict(case)
        for key in ("nod", "nlag", "lag", "azm", "displayed_dirs", "fixed_positions"):
            row[key] = join_ints(row[key])
        case_rows.append(row)
    case_fields = list(case_rows[0]) if case_rows else []
    write_csv(run_root / "case_matrix.csv", case_rows, case_fields)

    manifest = source_manifest(
        repo_root,
        [
            build_dir / "CMakeCache.txt",
            build_dir / "build.ninja",
            template_exe,
            window_exe,
            data_path,
        ]
        + ([paper_path] if paper_path is not None and paper_path.exists() else []),
    )
    write_csv(run_root / "sha256_manifest.csv", manifest, ["path", "size_bytes", "modified_local", "sha256"])

    raw_rows: list[dict] = []
    final_outputs: dict[tuple[str, str], Path] = {}
    total_routes = len(routes)
    for route_index, route_run in enumerate(routes, start=1):
        case = route_run.case
        target_property = f"{case.case_id}_{route_run.route}"
        case_output_dir = outputs_dir / target_property
        config_path = configs_dir / f"{target_property}.ini"
        config_text = (
            template_config(route_run, data_path, case_output_dir)
            if case.mode == "template"
            else window_config(route_run, data_path, case_output_dir)
        )
        config_path.write_text(config_text, encoding="utf-8")
        exe = template_exe if case.mode == "template" else window_exe
        phases = [("warmup", index + 1) for index in range(args.warmups)] + [
            ("recorded", index + 1) for index in range(args.repeats)
        ]
        print(
            f"[{route_index}/{total_routes}] {case.case_id} {route_run.route} "
            f"({len(phases)} executions)",
            flush=True,
        )
        for phase, run_number in phases:
            result = execute(exe, config_path, args.timeout)
            log_stem = f"{target_property}_{phase}_{run_number:02d}"
            (logs_dir / f"{log_stem}.stdout.txt").write_text(result.pop("stdout"), encoding="utf-8")
            (logs_dir / f"{log_stem}.stderr.txt").write_text(result.pop("stderr"), encoding="utf-8")
            status = (
                "TIMEOUT"
                if result["timed_out"]
                else ("PASS" if result["return_code"] == 0 else "FAILED")
            )
            raw_rows.append(
                {
                    "case_id": case.case_id,
                    "group": case.group,
                    "mode": case.mode,
                    "order_r": case.order_r,
                    "direction_count": case.direction_count,
                    "occurrence_vector": case.occurrence_vector,
                    "workload": case.workload,
                    "workload_label": case.workload_label,
                    "route": route_run.route,
                    "phase": phase,
                    "run_number": run_number,
                    "status": status,
                    "config": str(config_path),
                    **result,
                }
            )
            if status != "PASS":
                print(f"  {status}: see {logs_dir / (log_stem + '.stderr.txt')}", flush=True)
                break
        try:
            final_outputs[(case.case_id, route_run.route)] = output_file(case_output_dir, target_property)
        except FileNotFoundError:
            pass

    raw_fields = [
        "case_id",
        "group",
        "mode",
        "order_r",
        "direction_count",
        "occurrence_vector",
        "workload",
        "workload_label",
        "route",
        "phase",
        "run_number",
        "status",
        "return_code",
        "timed_out",
        "wall_time_s",
        "internal_time_s",
        "peak_working_set_mb",
        "peak_private_mb",
        "config",
    ]
    write_csv(run_root / "raw_runs.csv", raw_rows, raw_fields)

    comparisons: list[dict] = []
    for case in cases:
        bell_path = final_outputs.get((case.case_id, "bell"))
        avl_path = final_outputs.get((case.case_id, "avl"))
        if bell_path is None or avl_path is None:
            comparisons.append(
                {
                    "case_id": case.case_id,
                    "bell_count": 0,
                    "avl_count": 0,
                    "shared_count": 0,
                    "support_match": False,
                    "max_abs_error": math.inf,
                    "mean_abs_error": math.inf,
                    "rmse": math.inf,
                    "relative_l2": math.inf,
                    "nonfinite_value_count": -1,
                    "tolerance_failures": -1,
                    "pass": False,
                }
            )
            continue
        comparisons.append({"case_id": case.case_id, **compare_outputs(bell_path, avl_path)})

    comparison_fields = [
        "case_id",
        "bell_count",
        "avl_count",
        "shared_count",
        "support_match",
        "max_abs_error",
        "mean_abs_error",
        "rmse",
        "relative_l2",
        "nonfinite_value_count",
        "tolerance_failures",
        "pass",
    ]
    write_csv(run_root / "bell_avl_accuracy.csv", comparisons, comparison_fields)

    summaries = aggregate(raw_rows)
    summary_fields = [
        "case_id",
        "group",
        "mode",
        "order_r",
        "direction_count",
        "occurrence_vector",
        "workload",
        "workload_label",
        "route",
        "recorded_runs",
        "wall_time_s_mean",
        "wall_time_s_std",
        "wall_time_s_median",
        "internal_time_s_mean",
        "peak_working_set_mb_max",
        "peak_private_mb_max",
    ]
    write_csv(run_root / "route_summary.csv", summaries, summary_fields)

    pairs = build_pair_summary(summaries, comparisons)
    pair_fields = [
        "case_id",
        "group",
        "mode",
        "order_r",
        "direction_count",
        "occurrence_vector",
        "workload",
        "workload_label",
        "bell_wall_time_s_mean",
        "avl_wall_time_s_mean",
        "speedup_bell_over_avl",
        "bell_peak_working_set_mb",
        "avl_peak_working_set_mb",
        "memory_ratio_bell_over_avl",
        "bell_count",
        "avl_count",
        "shared_count",
        "support_match",
        "max_abs_error",
        "mean_abs_error",
        "rmse",
        "relative_l2",
        "nonfinite_value_count",
        "tolerance_failures",
        "pass",
    ]
    write_csv(run_root / "pair_summary.csv", pairs, pair_fields)

    all_processes_pass = all(row["status"] == "PASS" for row in raw_rows)
    all_accuracy_pass = all(bool(row["pass"]) for row in comparisons)
    gate = {
        "all_processes_pass": all_processes_pass,
        "all_bell_avl_accuracy_checks_pass": all_accuracy_pass,
        "overall_pass": all_processes_pass and all_accuracy_pass,
        "case_count": len(cases),
        "route_count": len(routes),
        "recorded_run_count": sum(row["phase"] == "recorded" for row in raw_rows),
        "failed_cases": [row["case_id"] for row in comparisons if not row["pass"]],
        "evidence_root": str(run_root),
    }
    (run_root / "gate.json").write_text(
        json.dumps(gate, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_run_report(run_root, pairs, gate)
    print(json.dumps(gate, indent=2, ensure_ascii=False), flush=True)
    return 0 if gate["overall_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
