from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from directional_validation import (
    connectivity_metrics,
    directional_offset,
    empirical_variogram,
    generate_accepted_iaaft_surrogates,
    generate_continuous_h2_surrogates,
    generate_discrete_h2_surrogates,
    h2_match_diagnostics,
    h2_offsets_from_templates,
    iaaft_surrogate,
    read_columnar_grid,
    rotate_offsets_90,
    rotation_exchange_error,
    signed_integer_direction,
    symmetric_leave_one_out_test,
    template_response,
    template_response_nd,
    transform_offsets_2d,
    validate_iaaft_surrogate,
    write_columnar_grid,
)


def test_rotation_of_explicit_offsets_and_direction_labels():
    offsets = np.array([[0, 0], [1, 0], [2, 0], [3, 0]])
    assert np.array_equal(rotate_offsets_90(offsets), np.array([[0, 0], [0, 1], [0, 2], [0, 3]]))
    assert directional_offset(90, 3) == (0, 3)


def test_iaaft_preserves_histogram_and_template_response_is_finite():
    field = np.arange(64, dtype=float).reshape(8, 8)
    surrogate = iaaft_surrogate(field, 4, 10)
    assert np.array_equal(np.sort(field.reshape(-1)), np.sort(surrogate.reshape(-1)))
    assert np.isfinite(template_response(field, np.array([[0, 0], [1, 0], [2, 0], [3, 0]])))
    assert rotation_exchange_error({0: 1.0, 45: 2.0}, {90: 1.0, 135: 2.0}) == {"0_to_90": 0.0, "45_to_135": 0.0}


def test_explicit_templates_support_integer_rotation_and_three_dimensional_offsets():
    line = np.array([[0, 0], [1, 0], [2, 0], [3, 0]])
    assert np.array_equal(transform_offsets_2d(line, 90, 2), np.array([[0, 0], [0, 2], [0, 4], [0, 6]]))
    volume = np.arange(5 * 6 * 7, dtype=float).reshape(5, 6, 7)
    tripod = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert np.isfinite(template_response_nd(volume, tripod))


def test_iaaft_acceptance_gate_uses_only_histogram_and_variogram():
    field = np.arange(64, dtype=float).reshape(8, 8)
    surrogate = iaaft_surrogate(field, 12, 20)
    verdict = validate_iaaft_surrogate(field, surrogate, [1, 2, 3], 100.0)
    assert verdict["histogram_exact"]
    assert verdict["accepted"]
    accepted, audit = generate_accepted_iaaft_surrogates(field, 2, 20, 15, [1, 2], 100.0, 2)
    assert len(accepted) == 2
    assert all(row["accepted"] for row in audit)
    assert np.isfinite(empirical_variogram(field, [1, 2])).all()


def test_connectivity_reports_largest_component_and_directional_runs():
    mask = np.array(
        [[True, True, False, False], [False, True, False, True], [False, False, False, True]],
        dtype=bool,
    )
    metrics = connectivity_metrics(mask, ("y", "x"))
    assert metrics["active_voxels"] == 5
    assert metrics["component_count"] == 2
    assert metrics["largest_component_voxels"] == 3
    assert metrics["run_max_x"] == 2
    assert metrics["run_max_y"] == 2


def test_columnar_grid_round_trip_preserves_three_dimensional_layout(tmp_path):
    values = np.arange(24, dtype=float).reshape(2, 3, 4)
    path = tmp_path / "volume.grid"
    write_columnar_grid(path, values, "porosity", x=[4, 5, 6, 7], y=[8, 9, 10], z=[11, 12])
    loaded = read_columnar_grid(path, "porosity")
    assert np.array_equal(loaded.values, values)
    assert np.array_equal(loaded.x, [4, 5, 6, 7])
    assert np.array_equal(loaded.y, [8, 9, 10])
    assert np.array_equal(loaded.z, [11, 12])


def test_signed_directions_cover_all_eight_orientations_without_modulo_reduction():
    expected = {
        0: (1, 0), 45: (1, 1), 90: (0, 1), 135: (-1, 1),
        180: (-1, 0), 225: (-1, -1), 270: (0, -1), 315: (1, -1),
    }
    assert {angle: signed_integer_direction(angle) for angle in expected} == expected
    assert signed_integer_direction(0) != signed_integer_direction(180)


def test_h2_offsets_include_every_pair_difference_and_its_reverse():
    templates = [
        np.array([[0, 0], [0, 1], [1, 0], [1, 0]], dtype=int),
        np.array([[0, 0], [2, 2], [2, 2], [0, 2]], dtype=int),
    ]
    offsets = {tuple(value) for value in h2_offsets_from_templates(templates)}
    assert (1, 0) in offsets
    assert (-1, 0) in offsets
    assert (2, 2) in offsets
    assert (-2, -2) in offsets
    assert (0, 0) in offsets


def test_h2_gate_and_surrogate_generators_preserve_their_frozen_inputs():
    field = np.arange(36, dtype=float).reshape(6, 6)
    offsets = np.array([[0, 0], [1, 0]], dtype=int)
    matching = h2_match_diagnostics(field, field.copy(), offsets)
    assert matching["accepted"]
    altered = field.copy()
    altered[0, 0] += 1000.0
    assert not h2_match_diagnostics(field, altered, [(0, 0), (1, 0)])["accepted"]

    continuous, continuous_audit = generate_continuous_h2_surrogates(
        field, 1, 8, 8, np.array([[0, 0]], dtype=int), 2, 0
    )
    assert len(continuous) == 1
    assert continuous_audit[-1]["accepted"]
    assert np.array_equal(np.sort(continuous[0].reshape(-1)), np.sort(field.reshape(-1)))

    occupancy = np.array(
        [[0.0, 0.125, 0.25], [0.375, 0.625, 0.75], [0.875, 1.0, 0.0]], dtype=float
    )
    discrete, discrete_audit = generate_discrete_h2_surrogates(
        occupancy, np.ones_like(occupancy, dtype=bool), 1, 9, offsets, 2, 20
    )
    assert len(discrete) == 1
    assert discrete_audit[-1]["accepted"]
    assert np.array_equal(np.sort(discrete[0].reshape(-1)), np.sort(occupancy.reshape(-1)))
    assert discrete_audit[-1]["difference_fraction_from_original"] > 0.0
    assert discrete_audit[-1]["occupancy_h2_rows"]
    assert discrete_audit[-1]["label_h2_rows"]
    with pytest.raises(ValueError, match="swap_iterations > 0"):
        generate_discrete_h2_surrogates(
            occupancy, np.ones_like(occupancy, dtype=bool), 1, 9, offsets, 2, 0
        )


def test_symmetric_surrogate_statistic_is_invariant_to_candidate_order():
    responses = np.array(
        [[1.0, 2.0, 4.0], [0.9, 2.2, 3.8], [1.2, 1.8, 4.1], [0.8, 2.1, 4.2]], dtype=float
    )
    baseline = symmetric_leave_one_out_test(responses)
    permutation = np.array([2, 0, 3, 1])
    permuted = symmetric_leave_one_out_test(responses[permutation])
    restored_scores = np.empty_like(permuted["scores"])
    restored_scores[permutation] = permuted["scores"]
    assert np.allclose(baseline["scores"], restored_scores)
    assert baseline["numerical_floor"] == permuted["numerical_floor"]
    original_index = int(np.flatnonzero(permutation == 0)[0])
    assert np.isclose(baseline["original_score"], permuted["scores"][original_index])
