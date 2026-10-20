from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
R1_5 = ROOT.parent / "r1_5_incremental_value_20260924"
for path in (ROOT, R1_5):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from prepare_stanford_vie import (
    PROPERTY_NAMES,
    build_manifest,
    coarse_volumes,
    validate_co_location_fixtures,
    write_subvolume,
)


def _synthetic_config(tmp_path: Path) -> dict[str, object]:
    porosity = np.arange(64, dtype=float)
    facies = np.zeros((4, 4, 4), dtype=int)
    facies[0, 0, :2] = 2
    facies[0, 1, :2] = 2
    facies[1, 0, :2] = 2
    # This separate fine 2-by-2-by-2 block has four channel cells and must
    # remain an invalid strict-majority tie after the common reversal.
    facies[2:4, 2:4, 0] = 2
    porosity_path = tmp_path / "poro.txt"
    facies_path = tmp_path / "facies.dat"
    porosity_path.write_text("h\n1\np\n" + "\n".join(map(str, porosity)), encoding="utf-8")
    facies_path.write_text("h\n1\nf\n" + " ".join(map(str, facies.reshape(-1))), encoding="utf-8")
    return {
        "fine_shape_zyx": [4, 4, 4],
        "porosity_path": str(porosity_path),
        "facies_path": str(facies_path),
        "porosity_header_lines": 3,
        "facies_header_lines": 3,
        "channel_code": 2,
    }


def test_coarse_volumes_reverse_axes_and_emit_colocated_occupancy(tmp_path: Path) -> None:
    config = _synthetic_config(tmp_path)
    volumes = coarse_volumes(config)
    porosity = np.arange(64, dtype=float)
    expected = porosity.reshape(4, 4, 4).reshape(2, 2, 2, 2, 2, 2).mean((1, 3, 5))[::-1, ::-1, :]
    assert np.array_equal(volumes.porosity, expected)
    assert volumes.channel_count.max() == 6
    assert np.array_equal(volumes.channel_occupancy, volumes.channel_count / 8.0)
    assert volumes.channel_strict_majority.sum() == 1
    assert np.count_nonzero(~volumes.channel_label_valid) == 1
    assert volumes.channel_count[0, 0, 0] == 4
    assert not volumes.channel_label_valid[0, 0, 0]


def test_fixed_coordinate_fixture_checks_all_properties_after_axis_reversal(tmp_path: Path) -> None:
    volumes = coarse_volumes(_synthetic_config(tmp_path))
    fixtures = [
        {
            "coordinate": [1, 2, 2],
            "expected": {
                "porosity": 10.5,
                "channel_count": 6,
                "channel_occupancy": 0.75,
                "channel_strict_majority": 1,
                "channel_label_valid": 1,
            },
        },
        {
            "coordinate": [1, 1, 1],
            "expected": {
                "porosity": 50.5,
                "channel_count": 4,
                "channel_occupancy": 0.5,
                "channel_strict_majority": 0,
                "channel_label_valid": 0,
            },
        },
    ]
    result = validate_co_location_fixtures(volumes, fixtures, porosity_atol=1.0e-15)
    assert len(result) == 2
    assert all(row["status"] == "PASS" for row in result)

    broken = [dict(fixtures[0], expected={**fixtures[0]["expected"], "channel_count": 0})]
    with pytest.raises(RuntimeError, match="channel_count mismatch"):
        validate_co_location_fixtures(volumes, broken, porosity_atol=1.0e-15)


def test_subvolume_manifest_records_each_property_row_count_and_hash(tmp_path: Path) -> None:
    config = _synthetic_config(tmp_path)
    volumes = coarse_volumes(config)
    _, subvolume = write_subvolume(tmp_path / "output", volumes, 1, 2)
    assert subvolume["shape_zyx"] == [2, 2, 2]
    assert subvolume["row_count"] == 8
    assert subvolume["tie_count"] == 1
    assert {row["property_name"] for row in subvolume["outputs"]} == set(PROPERTY_NAMES)
    assert all(row["row_count"] == 8 and len(row["sha256"]) == 64 for row in subvolume["outputs"])

    poro_inc_path = tmp_path / "poro.inc"
    poro_inc_path.write_text("PORO\n/\n", encoding="utf-8")
    config_path = tmp_path / "config.json"
    config.update(
        {
            "poro_inc_path": str(poro_inc_path),
            "coarse_spacing_xyz_m": [50.0, 50.0, 2.0],
        }
    )
    config_path.write_text(json.dumps(config), encoding="utf-8")
    fixture_path = tmp_path / "fixture.json"
    fixture_path.write_text(json.dumps({"fixtures": []}), encoding="utf-8")
    manifest = build_manifest(
        config=config,
        config_path=config_path,
        fixture_path=fixture_path,
        fixture_results=[],
        volumes=volumes,
        poro_inc_error=0.0,
        subvolume_metadata=[subvolume],
    )
    assert manifest["axis_mapping"]["input_array_order"] == ["z", "y", "x"]
    assert manifest["coarse_spacing_xyz_m"] == [50.0, 50.0, 2.0]
    assert manifest["subvolumes"][0]["outputs"][2]["property_name"] == "channel_occupancy"
