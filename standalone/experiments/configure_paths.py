"""Prepare machine-local direct-protocol configs without changing tracked files."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path


EXPERIMENTS = Path(__file__).resolve().parent
ROOT = EXPERIMENTS.parent


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve(source: Path, value: str) -> str:
    path = Path(value)
    return str((path if path.is_absolute() else source.parent / path).resolve())


def executable(build: Path, name: str) -> Path:
    filename = name + (".exe" if os.name == "nt" else "")
    for directory in (build / "Release", build / "bin", build):
        candidate = directory / filename
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"Build {name} in Release mode first: {build}")


def tool(value: str) -> str:
    found = shutil.which(value)
    if found:
        return str(Path(found).resolve())
    candidate = Path(value).resolve()
    if candidate.is_file():
        return str(candidate)
    raise FileNotFoundError(f"Tool not found: {value}; provide its full path or add it to PATH")


def configure(stanford_root: Path, build: Path, output: Path, cmake: str, ctest: str) -> Path:
    build, output = build.resolve(), output.resolve()
    inputs = {
        "porosity_path": stanford_root / "Porosity" / "porosity.txt",
        "facies_path": stanford_root / "Facies" / "facies.dat",
        "poro_inc_path": stanford_root / "Reservoir Simulation (Eclipse)" / "poro.inc",
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(f"Stanford VI-E source input is missing: {path}")
    if not (build / "CMakeCache.txt").is_file():
        raise FileNotFoundError(f"Configure the CMake build first: {build}")
    native = {name: executable(build, name) for name in (
        "template_cli", "movewindow_cli", "weighted_core_checks", "window_avl_core_checks"
    )}
    common_source = EXPERIMENTS / "configs" / "direct_common_protocol.json"
    common = load(common_source)
    common.update({
        "cmake_build_dir": str(build),
        "cmake_executable": tool(cmake),
        "ctest_executable": tool(ctest),
        "movewindow_cli": str(native["movewindow_cli"]),
        "template_cli": str(native["template_cli"]),
        "preprocessing_config": "preprocessing_config.json",
        "template_component_config": "template_component_config.json",
        "controlled_config": "controlled_config.json",
    })
    common["python_test_paths"] = [resolve(common_source, value) for value in common["python_test_paths"]]
    for dataset in common["datasets"]:
        for key in ("porosity", "label", "label_valid", "occupancy"):
            dataset[key] = resolve(common_source, dataset[key])

    preprocessing_source = EXPERIMENTS / "stanford_vie_preprocessing" / "preprocessing_config.json"
    preprocessing = load(preprocessing_source)
    preprocessing.update({key: str(path.resolve()) for key, path in inputs.items()})
    preprocessing["output_root"] = str((ROOT / "data" / "r1_5_high_order_20261010").resolve())
    preprocessing["co_location_fixture_path"] = resolve(preprocessing_source, preprocessing["co_location_fixture_path"])

    template_source = EXPERIMENTS / "template_directional_20261010" / "experiment_config.json"
    template = load(template_source)
    template["template_cli"] = str(native["template_cli"])
    template["native_test_workdir"] = str(build)
    template["native_test_commands"] = [
        {"name": "weighted_core_checks", "command": [str(native["weighted_core_checks"])]},
        {"name": "window_avl_core_checks", "command": [str(native["window_avl_core_checks"])]},
        {"name": "movewindow_cli_smoke", "command": [str(native["movewindow_cli"]), str(build / "tests" / "movewindow_cli_smoke.ini")]},
        {"name": "template_cli_exhaustive_integer_smoke", "command": [str(native["template_cli"]), str(build / "tests" / "template_exhaustive_smoke.ini")]},
    ]
    for key in ("output_root", "preprocessing_manifest"):
        template[key] = resolve(template_source, template[key])
    for dataset in template["datasets"]:
        for attribute in dataset["attributes"]:
            for key in ("grid_path", "valid_mask_path"):
                if key in attribute:
                    attribute[key] = resolve(template_source, attribute[key])

    payloads = {
        "direct_common_protocol.json": common,
        "preprocessing_config.json": preprocessing,
        "template_component_config.json": template,
        "controlled_config.json": load(EXPERIMENTS / "controlled_moment_experiments_20261010" / "experiment_config.json"),
    }
    for name in ("direct_section45.json", "direct_section46.json", "formal_section44.json"):
        payloads[name] = load(EXPERIMENTS / "configs" / name)
    # A fresh directory keeps previous run settings available for comparison.
    output.mkdir(parents=True, exist_ok=False)
    for name, payload in payloads.items():
        (output / name).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stanford-root", required=True, type=Path, help="Directory containing Porosity, Facies and Reservoir Simulation (Eclipse).")
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "local-config")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--ctest", default="ctest")
    args = parser.parse_args()
    result = configure(args.stanford_root, args.build_dir, args.output, args.cmake, args.ctest)
    print(f"Local configuration: {result}")
    print("No experiment has been started.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
