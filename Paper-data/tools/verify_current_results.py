"""Verify published result hashes, scopes, and figure/table provenance links."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=Path(__file__).resolve().parents[1] / "current_results/2026-10-12")
    package = parser.parse_args().package.resolve()
    copied = rows(package / "copied_files_manifest.csv")
    require(len({row["published_path"] for row in copied}) == len(copied), "Duplicate copy-manifest paths")
    for row in copied:
        path = (package / row["published_path"]).resolve()
        require(path.is_relative_to(package), f"Path escapes package: {path}")
        require(path.stat().st_size == int(row["bytes"]), f"Size mismatch: {path}")
        require(digest(path) == row["sha256"], f"Hash mismatch: {path}")

    derived = read_json(package / "figure_inputs/derived_manifest.json")["files"]
    for row in derived:
        path = package / "figure_inputs" / row["file"]
        require(path.stat().st_size == row["bytes"] and digest(path) == row["sha256"], f"Derived hash/size mismatch: {path}")
        require(len(rows(path)) == row["selected_rows"], f"Derived row-count mismatch: {path}")

    benchmark = package / "section_4_3"
    require(read_json(benchmark / "environment.json")["smoke"] is False, "Benchmark is a smoke run")
    require(read_json(benchmark / "gate.json")["overall_pass"] is True, "Benchmark gate failed")
    require(len(rows(benchmark / "case_matrix.csv")) == 16, "Benchmark case count")
    timings = rows(benchmark / "raw_runs.csv")
    require(sum(row["phase"] == "recorded" for row in timings) == 160, "Benchmark recorded run count")
    require(sum(row["phase"] == "warmup" for row in timings) == 32, "Benchmark warm-up count")
    direct = package / "sections_4_4_to_4_6"
    gate = read_json(direct / "gate.json")
    require(gate["mode"] == "formal" and gate["status"] == "PASS", "Direct run mode/status")
    require(gate["evaluation_mode"] == "direct_statistics", "Direct protocol mismatch")
    require(len(rows(direct / "cpp_feature_map_manifest.csv")) == gate["map_counts"]["map_count"] == 102, "Direct map count")
    component = read_json(direct / "section_4_4/run_20261011_191125/gate.json")
    require(component["status"] == "DIAGNOSTIC" and all(component["gates"].values()), "Section 4.4 component status/gates")
    require(component["surrogate_count"] == 199, "Section 4.4 surrogate count")

    contract = rows(direct / "report_artifact_contract.csv")
    for row in contract:
        folder = "figures" if row["artifact"].endswith(".png") else "tables"
        require(digest(direct / folder / row["artifact"]) == row["artifact_sha256"], f"Report artifact mismatch: {row['artifact']}")
    sources = read_json(package / "conditional_shape_supplement/source_hashes.json")
    for key in ("formal_evidence_manifest", "formal_gate", "frozen_common_protocol", "frozen_section45_protocol", "preprocessing_manifest", "run_manifest"):
        require(digest(direct / sources[key]["run_relative_path"]) == sources[key]["sha256"], f"Supplement provenance mismatch: {key}")
    require(sources["direct_analysis_values"]["sha256"] == derived[0]["source_sha256"], "Supplement/direct data provenance mismatch")
    associations = rows(package / "conditional_shape_supplement/conditional_rank_associations.csv")
    require(all(row["bootstrap_replicates"] == "2000" and row["rank_transform_scope"] == "complete_subvolume_fixed" for row in associations), "Conditional rank-association protocol mismatch")
    print(f"PASS: {len(copied)} exact-copy hashes, {len(derived)} derived CSV hashes/counts, {len(contract)} report contracts; formal run scopes and conditional-supplement provenance agree.")


if __name__ == "__main__":
    main()
