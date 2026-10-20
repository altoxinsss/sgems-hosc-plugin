# Section 4.4 Directional Template Run

Status: **DIAGNOSTIC**

## Gates

- `template_sampling_is_exhaustive_integer_offsets`: `PASS`
- `frozen_section44_protocol_binding`: `PASS`
- `preprocessing_manifest_inputs_verified`: `PASS`
- `preprocessing_axis_mapping_and_fixture_verified`: `PASS`
- `native_core_checks`: `PASS`
- `surrogate_count_complete`: `PASS`
- `surrogate_h2_gate`: `PASS`
- `surrogate_diversity`: `PASS`
- `bell_avl_max_absolute_error`: `PASS`
- `bell_avl_relative_l2_error`: `PASS`
- `rotation_signed_direction_max_absolute_error`: `PASS`
- `rotation_signed_direction_relative_l2_error`: `PASS`
- `required_figures`: `PASS`

## Frozen Design

- Production templates use signed integer offsets and exhaustive valid translations.
- Original fields are audited with Bell and AVL on identical template supports; surrogate fields use AVL only.
- The H2 surrogate gate reads only histogram and semivariogram quantities, not cumulant or connectivity results.
- Requested accepted surrogates per two-dimensional field: `199`.

## Output

- `tables/template_2d_response_long.csv` and `tables/template_3d_response_long.csv` are the production response facts.
- `tables/H2_offsets.csv`, `tables/directional_variogram_long.csv`, `tables/surrogate_acceptance.csv`, and `tables/surrogate_h2_diagnostics.csv` record the null-model gate.
- `tables/bell_avl_accuracy.csv`, `tables/rotation_validation_long.csv`, and `tables/connectivity_metrics.csv` record audit results.
- `evidence_manifest.json` and `sha256_manifest.csv` provide file-level provenance.
