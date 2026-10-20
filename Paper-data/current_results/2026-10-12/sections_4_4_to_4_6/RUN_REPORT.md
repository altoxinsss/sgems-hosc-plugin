# Direct-Statistics Sections 4.4-4.6 Run Report

Mode: `formal`
Status: `PASS`

## Direct Protocol

- Sections 4.5 and 4.6 use direct map statistics only; no classifier, cross-validation, learned threshold, fold, or OOF artifact is generated.
- C++ feature maps: `102`
- Primary maps: `90`
- Sensitivity maps: `12`
- lambda4/excess-kurtosis maximum absolute error: `2.663e-10`

## Checks

- `release_build`: PASS
- `python_regression_tests`: PASS
- `native_release_ctest`: PASS
- `stanford_preprocessing`: PASS
- `section44`: PASS
- `controlled_positive_controls`: PASS
- `cpp_python_parity_gate`: PASS
- `direct_cpp_map_matrix`: PASS
- `direct_report_artifacts`: PASS
- `formal_time_budget_projection`: PASS

## Core Artifacts

- `cpp_feature_map_manifest.csv`
- `direct_feature_summary.csv`
- `direct_effects.csv`
- `single_map_advantages.csv`
- `weighting_effects.csv`
- `bootstrap_summary.csv`
- `support_summary.csv`
- `direct_analysis_values.csv`
- `common_support_manifest.csv`
- `parameter_sensitivity_contrasts.csv`
- `lambda4_independent_audit.csv`
- `weighting_mechanism_diagnostics.csv`
- `time_projection.json`
- `report_artifact_contract.csv`
- `RUN_REPORT.md`
- `RESULTS_INDEX.md`
