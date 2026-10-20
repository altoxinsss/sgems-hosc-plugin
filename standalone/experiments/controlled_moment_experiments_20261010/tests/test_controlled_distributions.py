from __future__ import annotations

from decimal import Decimal

import numpy as np

from controlled_distributions import HIGH_THRESHOLD, FiniteDiscreteDistribution, controlled_experiments


def test_fourth_order_pair_matches_lower_cumulants_and_differs_at_fourth_order() -> None:
    fourth, _, _ = controlled_experiments()
    check = fourth.analytic_checks()
    assert check["matching_pass"]
    assert check["transformed_matching_pass"]
    assert check["high_probability_match"]
    assert check["difference_pass"]
    assert check["high_precision_matching_pass"]
    assert check["high_precision_transformed_matching_pass"]
    assert check["affine_scaling_pass"]
    assert check["higher_order_cumulants_reported"]
    assert check["analytic_gate_pass"]
    assert fourth.baseline_model == "M3"
    assert fourth.augmented_model == "M4"
    assert np.isclose(fourth.class_a.standardized_cumulant(4), -2.0)
    assert np.isclose(fourth.class_b.standardized_cumulant(4), 0.25)
    assert fourth.class_a.high_probability() == 0.5
    assert fourth.class_b.high_probability() == 0.5


def test_fifth_order_pair_matches_through_fourth_and_scales_nonzero_first_cumulant() -> None:
    _, fifth, _ = controlled_experiments()
    check = fifth.analytic_checks()
    assert check["matching_pass"]
    assert check["transformed_matching_pass"]
    assert check["high_probability_match"]
    assert check["difference_pass"]
    assert check["high_precision_matching_pass"]
    assert check["affine_scaling_pass"]
    assert check["analytic_gate_pass"]
    assert fifth.baseline_model == "M4"
    assert fifth.augmented_model == "M5"
    assert np.isclose(fifth.class_a.cumulants()[5], 32.0)
    assert np.isclose(fifth.class_b.cumulants()[5], 12.0)
    assert np.isclose(fifth.class_a.high_probability(), 23.0 / 30.0)
    assert np.isclose(fifth.class_b.high_probability(), 23.0 / 30.0)
    row_a = next(row for row in fifth.moment_rows() if row["class"] == "A")
    assert np.isclose(row_a["kappa1_x"], 0.5)
    assert np.isclose(row_a["kappa1_z"], 0.275)
    assert np.isclose(row_a["kappa5_z"], 0.05**5 * 32.0)
    direct = fifth.class_a.high_precision_direct_affine_cumulants()
    scaled = fifth.class_a.high_precision_affine_cumulants()
    assert abs(direct[1] - scaled[1]) < Decimal("1e-45")
    assert abs(direct[5] - scaled[5]) < Decimal("1e-45")


def test_sixth_order_pair_matches_through_fifth_and_differs_at_sixth_order() -> None:
    _, _, sixth = controlled_experiments()
    check = sixth.analytic_checks()
    assert check["matching_pass"]
    assert check["transformed_matching_pass"]
    assert check["high_probability_match"]
    assert check["difference_pass"]
    assert check["high_precision_matching_pass"]
    assert check["high_precision_transformed_matching_pass"]
    assert check["affine_scaling_pass"]
    assert check["analytic_gate_pass"]
    assert sixth.baseline_model == "M5"
    assert sixth.augmented_model == "M6"
    assert np.isclose(sixth.class_a.standardized_cumulant(6), 5.5)
    assert np.isclose(sixth.class_b.standardized_cumulant(6), 20.0 / 3.0)


def test_arbitrary_finite_distribution_supports_zero_and_asymmetric_values() -> None:
    distribution = FiniteDiscreteDistribution("asymmetric", (-2.0, 0.0, 3.0), (0.25, 0.5, 0.25))
    assert np.isclose(distribution.raw_moment(0), 1.0)
    assert np.isclose(distribution.high_probability(), 0.75)
    sampled = distribution.sample_windows(np.random.default_rng(123), 10, 27)
    assert np.all(np.isfinite(sampled))
    assert np.any(sampled >= HIGH_THRESHOLD)
    assert np.any(sampled < HIGH_THRESHOLD)
