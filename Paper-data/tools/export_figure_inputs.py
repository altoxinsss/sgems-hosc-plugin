"""Export small, lossless row selections from the archived formal run.

Uses only the Python standard library. Output CSV scalar text is preserved;
the exporter neither recomputes statistics nor changes the archived results.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export_rows(source, destination, keep, description, expected_hash):
    source_hash = sha256(source)
    if source_hash != expected_hash:
        raise ValueError(f"Source archive hash mismatch: {source}")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    input_count = output_count = 0
    with source.open(newline="", encoding="utf-8-sig") as input_file:
        reader = csv.DictReader(input_file)
        with destination.open("w", newline="", encoding="utf-8") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                input_count += 1
                if keep(row):
                    writer.writerow(row)
                    output_count += 1
    if not output_count:
        raise ValueError(f"Selection is empty: {source}")
    return {
        "file": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
        "source_file": str(source),
        "source_sha256": source_hash,
        "source_rows": input_count,
        "selected_rows": output_count,
        "selection": description,
        "transformation": "Row selection only; all columns and scalar strings preserved.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    source_run = args.source_run.resolve()
    output_dir = args.output_dir.resolve()
    with (source_run / "sha256_manifest.csv").open(newline="", encoding="utf-8-sig") as stream:
        expected = {row["path"].replace("\\", "/"): row["sha256"] for row in csv.DictReader(stream)}
    direct_relative = "direct_analysis_values.csv"
    template_relative = "section_4_4/run_20261011_191125/tables/template_2d_response_long.csv"
    records = [export_rows(
        source_run / direct_relative,
        output_dir / "direct_representative_z50.csv",
        lambda row: row["dataset_id"] == "z45_56"
        and row["window"] == "5x5x3" and row["z"] == "50"
        and row["parameter_role"] == "primary"
        and ((row["section"] == "4.5" and row["weighting"] == "uniform")
             or (row["section"] == "4.6" and row["weighting"] in
                 {"uniform", "coordinate_gaussian", "soft_high_porosity"})),
        "z45_56, z=50, 5x5x3, primary: Section 4.5 uniform and Section 4.6 three primary weightings",
        expected[direct_relative],
    ), export_rows(
        source_run / template_relative,
        output_dir / "template_2d_original_and_rotated.csv",
        lambda row: row["condition"] in {"original", "rotated90"},
        "All original and rotated90 rows, all datasets/properties/templates/algorithms; surrogate rows omitted",
        expected[template_relative],
    )]
    manifest = output_dir / "derived_manifest.json"
    if manifest.exists():
        raise FileExistsError(f"Refusing to overwrite: {manifest}")
    manifest.write_text(json.dumps({"schema_version": 1, "files": records}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
