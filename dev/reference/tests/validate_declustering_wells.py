"""Small integration audit on the existing 372 Stanford well samples.

This is a software check, not a new Section 4.5 effect experiment. The
kriging covariance below is a declared test setting, not a fitted model.
Outputs stay under the chosen build directory. Requires numpy only.
"""
import argparse
import configparser
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def run(executable, path):
    result = subprocess.run([str(executable), str(path)], capture_output=True, text=True, timeout=60)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    args = parser.parse_args()
    source_root = Path(__file__).resolve().parents[1]
    build = args.build_dir.resolve()
    executable = build / "Release/template_cli.exe"
    root = build / "test_output/stanford_three_methods"
    root.mkdir(parents=True, exist_ok=True)
    historical = source_root / "experiments/section44_45_stanford_20261110"
    source_ini = historical / "configs/wells_D_t2.5_z1_c15_o500.ini"
    settings = {
        "cell_declustering": None,
        "polygonal_declustering": {
            "dimensions": "2", "domain_min": "0.5,0.5", "domain_max": "75.5,100.5",
        },
        "kriging_declustering": {
            "dimensions": "2", "domain_min": "0.5,0.5", "domain_max": "75.5,100.5",
            "kriging_model": "exponential", "kriging_ranges": "15,15",
            "kriging_nugget_fraction": "0.1", "kriging_integration_cells": "75,100",
        },
    }
    method_results = {}
    support_reference = None
    for method, parameters in settings.items():
        for algorithm in ("AVL", "Bell"):
            config = configparser.ConfigParser(interpolation=None)
            config.read(source_ini, encoding="utf-8-sig")
            config["general"]["target_property"] = "result"
            config["general"]["calculate_algorithm"] = algorithm
            config["template_shape"]["nlag"] = "3,3,3,3"  # nine output cells only
            config["weighted_moment"]["weight_mode"] = method
            if parameters is not None:
                config["declustering"] = parameters
            output = root / (method + "_" + algorithm)
            config["runtime"]["output_dir"] = output.as_posix()
            config["runtime"]["write_template_tuples"] = "true"
            ini = root / (output.name + ".ini")
            with ini.open("w", encoding="utf-8") as stream:
                config.write(stream)
            run(executable, ini)
            audit = read_csv(output / "result.point_weights.csv")
            point_weights = np.array([float(row["weight"]) for row in audit])
            assert len(audit) == 372
            np.testing.assert_allclose(point_weights.sum(), 1, rtol=0, atol=2e-14)
            assert (point_weights > 0).all()
            support = read_csv(output / "result.support.csv")
            signature = [(r["lag_positions"], r["support_count"], r["anchor_count"]) for r in support]
            if support_reference is None:
                support_reference = signature
            assert signature == support_reference, "Weight selection changed complete-tuple support"
            tuples = read_csv(output / "result.tuples.csv")
            if method == "cell_declustering" and algorithm == "AVL":
                reference_path = historical / "results/wells_D_t2.5_z1_c15_o500/wells_D_t2.5_z1_c15_o500.support.csv"
                reference = {r["lag_positions"]: r for r in read_csv(reference_path)}
                for row in support:
                    old = reference[row["lag_positions"]]
                    for column in ("value", "weight_sum", "ESS"):
                        np.testing.assert_allclose(float(row[column]), float(old[column]), rtol=2e-13, atol=1e-18, equal_nan=True)
            values = np.array([float(row["value"]) for row in support])
            if algorithm == "AVL":
                method_results[method] = {"weights": point_weights, "values": values, "audit": audit}
            else:
                np.testing.assert_allclose(values, method_results[method]["values"], rtol=1e-10, atol=1e-16, equal_nan=True)
            # Independent raw joint moments + enumerated set partitions (order 5).
            for row in support:
                selected = [t for t in tuples if t["lag_cell"] == row["lag_cell"]]
                if not selected:
                    assert np.isnan(float(row["value"]))
                    continue
                x = np.array([[float(v) for v in t["values"].split(";")] for t in selected])
                weight = np.array([float(t["tuple_weight"]) for t in selected])
                for t, w in zip(selected, weight):
                    ids = sorted(set(int(v) - 1 for v in t["node_ids"].split(";")))
                    expected = np.exp(np.log(point_weights[ids]).mean())
                    np.testing.assert_allclose(w, expected, rtol=1e-13, atol=0)
                weight /= weight.sum()
                terms = []
                def partition(index, blocks):
                    if index == 5:
                        import math
                        factor = (-1) ** (len(blocks) - 1) * math.factorial(len(blocks) - 1)
                        terms.append(factor * np.prod([weight @ np.prod(x[:, block], axis=1) for block in blocks]))
                        return
                    for i in range(len(blocks)):
                        partition(index + 1, [b + [index] if j == i else b for j, b in enumerate(blocks)])
                    partition(index + 1, blocks + [[index]])
                partition(0, [])
                np.testing.assert_allclose(float(row["value"]), sum(terms), rtol=1e-9, atol=1e-15)

    audit = method_results["kriging_declustering"]["audit"]
    xy = np.array([[float(r["x"]), float(r["y"])] for r in audit])
    sites, inverse, counts = np.unique(xy, axis=0, return_inverse=True, return_counts=True)
    assert len(sites) == 31 and (counts == 12).all()
    # Independent augmented ordinary-kriging system, rather than Cholesky
    # elimination used by C++. Same explicit domain midpoint quadrature.
    xx, yy = np.meshgrid(np.arange(1, 76), np.arange(1, 101))
    domain = np.column_stack([xx.ravel(), yy.ravel()])
    distances = np.linalg.norm((sites[:, None] - sites[None, :]) / 15, axis=2)
    covariance = .9 * np.exp(-3 * distances) + .1 * np.eye(len(sites))
    b = (.9 * np.exp(-3 * np.linalg.norm((sites[:, None] - domain[None, :]) / 15, axis=2))).mean(axis=1)
    system = np.block([[covariance, np.ones((31, 1))], [np.ones((1, 31)), np.zeros((1, 1))]])
    solution = np.linalg.solve(system, np.r_[b, 1])
    np.testing.assert_allclose(method_results["kriging_declustering"]["weights"], solution[:-1][inverse]/counts[inverse], rtol=1e-11, atol=1e-14)

    # Independent fine nearest-site integration checks polygon areas. This is
    # an approximate cross-check, separate from the exact analytic C++ tests.
    pixel_counts = np.zeros(31, dtype=np.int64)
    x = .5 + (np.arange(1000)+.5)*75/1000
    for y in .5 + (np.arange(1000)+.5)*100/1000:
        d = (x[:, None]-sites[None, :, 0])**2 + (y-sites[None, :, 1])**2
        pixel_counts += np.bincount(d.argmin(axis=1), minlength=31)
    exact_areas = np.bincount(inverse, weights=method_results["polygonal_declustering"]["weights"])
    polygon_error = float(np.max(np.abs(exact_areas-pixel_counts/1e6)))
    assert polygon_error < .0002, polygon_error
    report = {
        "purpose": "implementation validation only; covariance is not fitted manuscript evidence",
        "source_rows": 372, "distinct_xy_wells": 31, "lag_cells_per_run": 9,
        "valid_lag_cells": int(np.isfinite(method_results["cell_declustering"]["values"]).sum()),
        "same_support_all_methods": True, "historical_cell_regression": "passed",
        "bell_avl_and_independent_raw_partition": "passed",
        "kriging_independent_augmented_solve": "passed",
        "polygon_max_area_fraction_difference_to_1m_pixel_reference": polygon_error,
        "template_cli_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "methods": {m: {"min_point_weight": float(r["weights"].min()), "max_point_weight": float(r["weights"].max())}
                    for m, r in method_results.items()},
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
