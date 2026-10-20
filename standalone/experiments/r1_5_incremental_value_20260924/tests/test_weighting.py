import numpy as np

from statistics import compute_local_statistics
from weighting import (
    WeightConfig,
    coordinate_gaussian_window_weights,
    compute_weighted_statistics,
    generate_weights,
    joint_configuration_weights,
    matched_ess_uniform_weights,
    value_gaussian_weights,
)
from window_support import WindowBatch


def make_batch(values):
    values = np.asarray(values, dtype=float)
    return WindowBatch(values, np.array([0]), np.array([0]), np.array([0]), np.zeros_like(values, dtype=int))


def test_gaussian_weight_formula_and_minimum():
    values = np.array([[0.26, 0.34, 0.50]])
    config = WeightConfig(reference_value=0.26, bandwidth=0.08, min_weight=0.1)
    weights = value_gaussian_weights(values, config)
    assert np.isclose(weights[0, 0], 1.0)
    assert np.isclose(weights[0, 1], np.exp(-0.5))
    assert weights[0, 2] == 0.0


def test_uniform_weighting_matches_unweighted():
    windows = make_batch([[0.1, 0.2, 0.4, 0.5]])
    unweighted = compute_local_statistics(windows, 0.25, 1e-15)
    weighted = compute_weighted_statistics(windows, np.ones_like(windows.values), 0.25, 1e-15).statistics
    for name in unweighted.as_dict():
        assert np.allclose(unweighted.as_dict()[name], weighted.as_dict()[name], equal_nan=True)


def test_coordinate_gaussian_center_and_relative_window_distance():
    weights = coordinate_gaussian_window_weights((3, 3, 3), (2.0, 3.0, 4.0), bandwidth=2.0).reshape(3, 3, 3)
    assert np.isclose(weights[1, 1, 1], 1.0)
    # For a 3-cell axis, an edge is one normalized half-width from centre;
    # sigma_s=2 therefore gives exp(-1/(2*2^2)).
    assert np.isclose(weights[1, 1, 0], np.exp(-0.125))


def test_coordinate_gaussian_is_window_position_invariant():
    config = WeightConfig(mode="coordinate_gaussian", bandwidth=2.0, min_weight=0.0)
    first = generate_weights(np.arange(27, dtype=float).reshape(1, 27), config, window_size=(3, 3, 3))
    translated = generate_weights((1000.0 + np.arange(27)).reshape(1, 27), config, window_size=(3, 3, 3))
    assert np.allclose(first, translated)


def test_zero_weight_window_is_nan():
    windows = make_batch([[1, 2, 3]])
    result = compute_weighted_statistics(windows, np.zeros_like(windows.values), 0.25, 1e-15)
    assert np.isnan(result.statistics.mean[0])
    assert np.isnan(result.effective_sample_size[0])


def test_explicit_scheme_switches_are_mutually_exclusive():
    config = WeightConfig(enable_fixed_value_gaussian=True, enable_soft_high_porosity=True)
    try:
        config.validate()
    except ValueError as error:
        assert "Only one" in str(error)
    else:
        raise AssertionError("Conflicting switches must be rejected")


def test_soft_high_porosity_is_monotone_and_half_at_threshold():
    config = WeightConfig(
        enable_fixed_value_gaussian=False,
        enable_soft_high_porosity=True,
        high_porosity_threshold=0.25,
        soft_transition_width=0.02,
        soft_weight_strength=1.0,
    )
    weights = generate_weights(np.array([[0.20, 0.25, 0.30]]), config)
    assert np.isclose(weights[0, 1], 0.5)
    assert np.all(np.diff(weights[0]) > 0.0)


def test_adaptive_upper_tail_uses_each_window_context():
    config = WeightConfig(
        enable_fixed_value_gaussian=False,
        enable_adaptive_upper_tail=True,
        adaptive_global_threshold=0.0,
        adaptive_mean_offset=0.0,
        adaptive_transition_scale=1.0,
        adaptive_weight_strength=1.0,
    )
    values = np.array([[0.0, 1.0, 2.0], [10.0, 11.0, 12.0]])
    weights = generate_weights(values, config)
    assert np.allclose(weights[:, 1], 0.5)
    assert np.all(np.diff(weights, axis=1) > 0.0)


def test_joint_configuration_weight_is_geometric_mean():
    component_weights = np.array([[1.0, np.exp(-0.5)], [0.25, 1.0]])
    result = joint_configuration_weights(component_weights)
    assert np.allclose(result, [np.exp(-0.25), 0.5])


def test_scheme_identifier_is_traceable():
    config = WeightConfig(
        enable_fixed_value_gaussian=True,
        enable_joint_geometric_mean=True,
    )
    assert config.scheme_id() == "fixed_value_gaussian.v1+joint_geometric_mean.v1"


def test_formal_and_legacy_weighting_schemes_are_explicitly_labelled():
    assert WeightConfig(mode="uniform").scheme_status() == "formal"
    assert WeightConfig(mode="coordinate_gaussian").scheme_status() == "formal"
    assert WeightConfig(mode="soft_high_porosity").scheme_status() == "formal"
    assert WeightConfig(mode="value_gaussian").scheme_status() == "legacy"
    assert WeightConfig(mode="adaptive_upper_tail").scheme_status() == "legacy"
    assert WeightConfig(enable_fixed_value_gaussian=True, enable_joint_geometric_mean=True).scheme_status() == "legacy"


def test_matched_ess_uniform_control_is_reproducible_and_matches_expected_ess():
    reference = np.array([[1.0, 1.0, 0.0, 0.0], [1.0, 0.5, 0.5, 0.0]])
    one = matched_ess_uniform_weights(reference, 9)
    two = matched_ess_uniform_weights(reference, 9)
    assert np.array_equal(one, two)
    assert one[0].sum() == 2.0
    assert one[1].sum() in {2.0, 3.0}
