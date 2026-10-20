from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from run_template_directional_experiment import (
    _integer_direction,
    _matching_explicit_surface,
    _module_only_status,
    _nominal_physical_lag,
    _relative_l2,
    _write_cli_ini,
    _write_module_results_index,
)


def test_cli_ini_captures_production_directional_parameters(tmp_path):
    ini_path = tmp_path / "case.ini"
    _write_cli_ini(
        ini_path,
        tmp_path / "source.grid",
        "porosity",
        "fourth_order_case",
        tmp_path / "output",
        nod=[2, 1],
        nlag=[4, 4],
        lag=[1, 1],
        azimuth=[45, 135],
        dip=[0, 0],
        displayed_dirs=[1, 2],
        coordinate_dimensions=2,
        source_spacing=[50.0, 50.0, 2.0],
        direction_vectors=[[1, 1], [-1, 1]],
    )
    text = ini_path.read_text(encoding="utf-8")
    assert "calculate_algorithm=AVL Based algorithm" in text
    assert "template_sampling=exhaustive_integer_offsets" in text
    assert "integer_direction_vectors=1,1;-1,1" in text
    assert "NOD=2,1" in text
    assert "AZM=45,135" in text
    assert "displayed_dirs=1,2" in text


def test_cli_ini_captures_exact_three_dimensional_xyz_geometry(tmp_path):
    ini_path = tmp_path / "xyz.ini"
    _write_cli_ini(
        ini_path,
        tmp_path / "source.grid",
        "porosity",
        "fourth_order_xyz",
        tmp_path / "output",
        nod=[1, 1, 1],
        nlag=[4, 4, 4],
        lag=[1, 1, 1],
        azimuth=[0, 90, 0],
        dip=[0, 0, 90],
        displayed_dirs=[1, 2, 3],
        coordinate_dimensions=3,
        source_spacing=[50.0, 50.0, 2.0],
        direction_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    text = ini_path.read_text(encoding="utf-8")
    assert "coordinate_dimensions=3" in text
    assert "integer_direction_vectors=1,0,0;0,1,0;0,0,1" in text
    assert "NOD=1,1,1" in text
    assert "NLAG=4,4,4" in text
    assert "displayed_dirs=1,2,3" in text


def test_module_only_path_cannot_emit_final_status_or_formal_result_index(tmp_path):
    assert _module_only_status({"all_checks": True}) == "DIAGNOSTIC"
    assert _module_only_status({"one_failed_check": False}) == "FAIL"
    index = _write_module_results_index(tmp_path)
    assert index.name == "MODULE_RESULTS_INDEX.md"
    assert not (tmp_path / "RESULTS_INDEX.md").exists()
    assert "not a formal Section 4.4 evidence package" in index.read_text(encoding="utf-8")


def test_rotation_error_and_nominal_lag_metadata_are_deterministic():
    import numpy as np

    assert _relative_l2(np.ones((2, 2)), np.ones((2, 2))) == 0.0
    assert _relative_l2(np.ones((2, 2)), np.ones((2, 3))) == float("inf")
    lag = _nominal_physical_lag(0, 0, 2, [3.0, 4.0, 5.0])
    assert np.isclose(lag["nominal_dx"], 6.0)
    assert np.isclose(lag["nominal_dy"], 0.0)
    assert np.isclose(lag["nominal_distance"], 6.0)


def test_matching_explicit_surface_uses_the_configured_nod_2_1_geometry():
    import numpy as np

    field = np.arange(64, dtype=float).reshape(8, 8)
    case = {"nod": [2, 1], "azimuth": [0, 90], "lag_cells": [1, 1]}
    surface = _matching_explicit_surface(field, case, np.zeros((3, 4)))
    assert surface.shape == (3, 4)
    assert np.isfinite(surface).all()
    assert _integer_direction(0) == (1, 0)
    assert _integer_direction(90) == (0, 1)
