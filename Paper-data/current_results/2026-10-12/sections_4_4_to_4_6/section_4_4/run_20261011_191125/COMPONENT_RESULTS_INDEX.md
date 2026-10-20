# Section 4.4 Component Result Placement

This component index is consumed by the top-level formal entrypoint. The study-level result index and final PASS status are written only at the top level.

## 4.4.1 Directional integer-offset responses

Use `tables/template_2d_response_long.csv`, `tables/template_3d_response_long.csv`, and `tables/template_definitions.csv`. The response rows retain every lag cell, signed integer offset, physical offset, and complete-translation support count.

## 4.4.2 Second-order-matched surrogate comparison

Use `tables/H2_offsets.csv`, `tables/directional_variogram_long.csv`, `tables/surrogate_acceptance.csv`, `tables/surrogate_h2_diagnostics.csv`, `tables/surrogate_response_summary.csv`, and `tables/multiplicity_decisions.csv`. The porosity omnibus is the primary response test; occupancy and strict-majority connectivity are separate Holm-adjusted secondary tests.

## 4.4.3 Rotation, algorithm audit, and connectivity

Use `tables/rotation_validation_long.csv`, `tables/bell_avl_accuracy.csv`, and `tables/connectivity_metrics.csv`. Explicit line/rectangle/L responses are marked exploratory and do not enter the omnibus surrogate tests.

## Reproducibility

Use generated `configs/`, `inputs/`, `logs/`, `sha256_manifest.csv`, and `evidence_manifest.json`.
