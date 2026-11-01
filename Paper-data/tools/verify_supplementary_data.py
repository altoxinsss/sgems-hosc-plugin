"""Verify the current Data S1 ZIP and its complete member checksum inventory."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "supplementary_data_s1",
    )
    package = parser.parse_args().package.resolve()
    manifest = json.loads((package / "archive_manifest.json").read_text(encoding="utf-8"))
    archive = package / manifest["filename"]
    payload = archive.read_bytes()
    require(len(payload) == manifest["bytes"], "Archive byte-count mismatch")
    require(hashlib.sha256(payload).hexdigest() == manifest["sha256"], "Archive hash mismatch")
    with zipfile.ZipFile(archive) as zipped:
        names = zipped.namelist()
        require(len(names) == len(set(names)), "Duplicate archive members")
        require(len(names) == manifest["members"], "Archive member-count mismatch")
        for name in names:
            path = PurePosixPath(name)
            require(not path.is_absolute() and ".." not in path.parts and "\\" not in name,
                    f"Unsafe member path: {name}")
        expected = {}
        for line in zipped.read("SHA256SUMS.txt").decode("utf-8-sig").splitlines():
            checksum, name = line.split("  ", 1)
            require(name not in expected, f"Duplicate checksum entry: {name}")
            expected[name] = checksum
        require(set(expected) == set(names) - {"SHA256SUMS.txt"},
                "Checksum inventory does not cover exactly the package members")
        for name, checksum in expected.items():
            actual = hashlib.sha256(zipped.read(name)).hexdigest()
            require(actual == checksum, f"Member hash mismatch: {name}")
        required = {
            "TABLE_TO_FILE_MAP.csv", "FIGURE_TO_FILE_MAP.csv",
            "matched_fourth_order_windows/scripts/verify_selected_windows.py",
            "matched_fourth_order_windows/scripts/rebuild_screening.py",
            "matched_fourth_order_windows/Table_S11_matched_attributes.csv",
            "matched_fourth_order_windows/figures/Figure_S16_strict_matched_windows.svg",
            "matched_fourth_order_windows/figures/Figure_S16_strict_matched_windows.tiff",
        }
        require(required <= set(names), "Required current manuscript evidence is missing")
    print(f"PASS: archive SHA-256, {len(names)} members, and {len(expected)} member checksums agree.")


if __name__ == "__main__":
    main()
