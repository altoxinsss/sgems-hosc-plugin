import numpy as np

from statistics import compute_local_statistics
from window_support import WindowBatch


def batch(values):
    values = np.asarray(values, dtype=float)
    n = len(values)
    return WindowBatch(values, np.arange(n), np.zeros(n), np.zeros(n), np.zeros_like(values, dtype=int))


def test_known_moments_and_cumulant():
    result = compute_local_statistics(batch([[1, 2, 3, 4]]), 3.0, 1e-15)
    assert np.allclose(result.mean, [2.5])
    assert np.allclose(result.variance, [1.25])
    assert np.allclose(result.central_moment4, [2.5625])
    assert np.allclose(result.cumulant4, [-2.125])
    assert np.allclose(result.excess_kurtosis, [-1.36])
    assert np.allclose(result.high_fraction, [0.5])


def test_zero_variance_kurtosis_is_nan():
    result = compute_local_statistics(batch([[2, 2, 2]]), 1.0, 1e-15)
    assert np.isnan(result.excess_kurtosis[0])
    assert result.cumulant4[0] == 0.0


def test_fifth_and_sixth_cumulant_explicit_identities():
    windows = batch([[0.1, 0.2, 0.4, 0.7, 1.1]])
    centered = windows.values - np.mean(windows.values, axis=1)[:, None]
    m2 = np.mean(centered**2, axis=1)
    m3 = np.mean(centered**3, axis=1)
    m4 = np.mean(centered**4, axis=1)
    expected5 = np.mean(centered**5, axis=1) - 10 * m3 * m2
    expected6 = np.mean(centered**6, axis=1) - 15 * m4 * m2 - 10 * m3**2 + 30 * m2**3
    result5 = compute_local_statistics(windows, 0.25, 1e-15, 5)
    result6 = compute_local_statistics(windows, 0.25, 1e-15, 6)
    assert np.allclose(result5.as_dict()["cumulant5"], expected5)
    assert np.allclose(result6.as_dict()["cumulant6"], expected6)


def test_complete_statistic_set_marks_low_variance_invalid_and_kurtosis_equivalence():
    result = compute_local_statistics(batch([[0.0, 1.0, 2.0, 3.0], [2.0, 2.0, 2.0, 2.0]]), 2.0, 1e-12)
    fields = result.as_dict()
    assert {"cumulant2", "cumulant3", "cumulant4", "cumulant5", "cumulant6", "skewness", "standardized_kappa4", "standardized_kappa5", "standardized_kappa6", "valid"}.issubset(fields)
    assert np.isclose(fields["excess_kurtosis"][0], fields["cumulant4"][0] / fields["variance"][0] ** 2)
    assert np.allclose(fields["standardized_kappa4"], fields["excess_kurtosis"], equal_nan=True)
    assert fields["valid"][0]
    assert not fields["valid"][1]
    assert np.isnan(fields["standardized_kappa6"][1])
