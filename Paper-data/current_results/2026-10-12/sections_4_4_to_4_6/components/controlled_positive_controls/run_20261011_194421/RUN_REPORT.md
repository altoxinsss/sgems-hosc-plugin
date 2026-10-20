# Direct Matching-Moment Control Report

Run type: `formal`.

The control compares finite-window cumulant estimators directly. It does not fit a classifier, split training and test samples, or compute predictive metrics.

## Target-order estimates

| Pair | Samples | Order | Analytic difference | Estimated difference | 95% CI | Detected |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| fifth_order_matching_moments | 27 | 5 | 6.25e-06 | 3.1652373e-06 | [3.1184212e-06, 3.2126409e-06] | True |
| fifth_order_matching_moments | 75 | 5 | 6.25e-06 | 5.0518237e-06 | [5.0202397e-06, 5.0848274e-06] | True |
| fifth_order_matching_moments | 147 | 5 | 6.25e-06 | 5.6308864e-06 | [5.6025433e-06, 5.6573938e-06] | True |
| fourth_order_matching_moments | 27 | 4 | -1.40625e-05 | -1.0644856e-05 | [-1.0679795e-05, -1.0609895e-05] | True |
| fourth_order_matching_moments | 75 | 4 | -1.40625e-05 | -1.2783416e-05 | [-1.2804665e-05, -1.2761128e-05] | True |
| fourth_order_matching_moments | 147 | 4 | -1.40625e-05 | -1.3399506e-05 | [-1.3409541e-05, -1.3389506e-05] | True |
| sixth_order_matching_moments | 27 | 6 | -1.8229167e-08 | -1.4651052e-09 | [-2.410724e-09, -5.0862592e-10] | True |
| sixth_order_matching_moments | 75 | 6 | -1.8229167e-08 | -1.1183284e-08 | [-1.1613198e-08, -1.0756286e-08] | True |
| sixth_order_matching_moments | 147 | 6 | -1.8229167e-08 | -1.4630997e-08 | [-1.4897061e-08, -1.4368419e-08] | True |

## Gates

| Gate | Value |
| --- | --- |
| analytic_moment_checks_pass | True |
| expected_result_rows | 1620 |
| actual_result_rows | 1620 |
| complete_execution_pass | True |
| all_metrics_finite | True |
| target_order_detection_pass | True |
| matched_lower_order_diagnostic_mean_abs_delta | 1.0187553860263997e-05 |
| overall_pass | True |
| status | VALIDATION_PASS |
| run_kind | formal |
| component_execution | COMPLETE |

## Files

- `population_moments.csv`: analytic cumulants for the paired distributions.
- `direct_per_seed_estimates.csv`: finite-window direct cumulant estimates by seed.
- `direct_summary.csv`: bias, RMSE, group differences, and confidence intervals.
- `gate.json`: execution and target-order detection gates.
