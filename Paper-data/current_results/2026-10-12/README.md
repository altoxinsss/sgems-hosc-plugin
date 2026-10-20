# October 2026 result snapshot

This package contains selected results supporting Sections 4.3-4.6 of *An integrated SGeMS plugin for computing high-order spatial statistics*. It was curated on 2026-10-20 from completed local runs. Numerical results were copied without recomputation. The figure-input CSVs preserve all scalar text and columns from selected source rows.

## Source runs and status

The source base is `D:/Development/DLL/DLL/WithoutSGeMS/output/` on the recording machine.

| Published directory | Source path relative to the source base | Recorded status |
| --- | --- | --- |
| `section_4_3/` | `r1_11_scalability_20260909/run_20261012_084358` | `smoke=false`, `overall_pass=true` |
| `sections_4_4_to_4_6/` | `direct_4_4_4_6_20261011/20261012_031025/run_20261011_191032` | `mode=formal`, `status=PASS` |
| `sections_4_4_to_4_6/section_4_4/run_20261011_191125/` | Section 4.4 component within the direct run | `mode=formal`, `status=DIAGNOSTIC`; all 13 component gates true |
| `sections_4_4_to_4_6/components/controlled_positive_controls/run_20261011_194421/` | Controlled component within the direct run | `run_kind=formal`, `status=VALIDATION_PASS` |
| `conditional_shape_supplement/` | `conditional_high_order_window_analysis_fixed_rank_20261012` | Completed post-processing of the formal direct run; provenance and bootstrap diagnostics in `RUN_REPORT.md`; no separate gate status |

The `20261012_031025` directory identifies the selected complete archive; its nested run ID remains `run_20261011_191032`. The direct entrypoint ID is `direct-44-46-374de37d-a8dc-49ad-a503-9ed66bc0c9ea-efe9b4ee7ea3`. Section 4.4 uses the component-only status `DIAGNOSTIC`; the complete study status is written by the top-level runner. Preserve both status levels when describing the evidence.

## Evidence and interpretation

| Scope | Principal files | Interpretation |
| --- | --- | --- |
| 4.3 performance and agreement | `section_4_3/{case_matrix,raw_runs,route_summary,pair_summary,bell_avl_accuracy}.csv` | 16 cases, 32 routes, one warm-up and five measured repeats per route; 160 recorded runs. Includes orders 4-8 and separate template/window cases. Memory is peak recorded process working set in MiB. Timing is specific to the recorded machine/build. |
| 4.4 directional templates | Component `tables/{template_definitions,template_3d_response_long,surrogate_response_summary,multiplicity_decisions}.csv` and `figure_inputs/template_2d_original_and_rotated.csv` | Exhaustive valid integer-offset translations. There are 199 accepted surrogates per 2D field. The porosity omnibus is primary; occupancy and strict-majority connectivity are Holm-adjusted secondary tests. Explicit line/rectangle/L templates are exploratory. |
| 4.4 audits and support | Component `tables/{bell_avl_accuracy,rotation_validation_long,connectivity_metrics,dataset_inventory,H2_offsets,surrogate_diversity}.csv` | Numerical agreement, rotation, supports, connectivity, and null-model design. Full H2 candidate/acceptance histories remain in the source archive. |
| 4.5 direct uniform statistics | `direct_effects.csv`, `single_map_advantages.csv`, `tables/section45_direct_effects_primary_window.csv` | Per-map contrasts and separation against simple baselines on the three slabs; primary window 5x5x3, scale windows 3x3x3 and 7x7x3. |
| 4.6 weighting | `weighting_effects.csv`, `parameter_sensitivity_contrasts.csv`, `tables/section46_*.csv` | Primary uniform, coordinate-Gaussian, and soft-high-porosity maps; paired weighting effects on common support. Parameter sensitivity is scoped to z45_56. |
| Uncertainty and validity | `bootstrap_summary.csv`, `support_summary.csv`, `common_support_manifest.csv`, `lambda4_independent_audit.csv` | Paired spatial-block bootstrap and exclusions. Preserve section, dataset, window, weighting, parameter set, feature, support, and bootstrap scope when selecting rows. |
| Positive controls and implementation checks | `components/controlled_positive_controls/`, `cpp_python_fixture_matrix/`, `cpp_python_parity/`, `uniform_three_route_audit/` | Analytic/finite-window controls and implementation agreement. These checks have their own scope and are distinct from Stanford-data effect estimates. |

Paths in the last four rows are relative to `sections_4_4_to_4_6/`. The current Sections 4.5/4.6 protocol is `direct_statistics`: classifier, cross-validation, OOF, and M4/M6 prediction metrics are not part of this run. No smoke run was selected for this package.

`effect_E` is the channel-minus-nonchannel mean difference divided by the square root of the mean of the two sample variances; `separation_S = abs(effect_E)`. `advantage_A` is feature separation minus the strongest applicable simple-baseline separation. `W` is weighted-minus-uniform separation; `G` is weighted-minus-uniform single-map advantage. `lambda_r = kappa_r / kappa2^(r/2)`; lambda4 is excess kurtosis. Negative advantages and confidence intervals crossing zero are retained and must be interpreted as reported.

The primary bootstrap uses 2,000 replicates and spatial blocks 8x13x6. Scale and parameter-sensitivity scopes use 1,000 replicates; the larger 12x20x6 block sensitivity uses 500. `macro_equal_dataset` denotes equal weighting of the three dataset estimates. The controlled component separately uses 30 seeds, 2,000 windows per class, and 10,000 bootstrap replicates.

## Conditional marginal and tail-shape supplement

`conditional_shape_supplement/` contains the later lower-order-controlled analysis used for the conditional shape comparison in Section 4.5. The recorded `source_hashes.json` binds it to the same formal direct table, porosity grids, protocol, and gate. It uses all 68,160 numerically valid uniform 5x5x3 windows per slab, independently of channel labels and facies-validity masks.

For lambda5 the rank projection controls for mean, log variance, high-porosity proportion, lambda3, and lambda4. For lambda6 it also controls for lambda5. `conditional_rank_associations.csv` reports correlations between the two residualized rank variables, confidence intervals, and slopes. The complete-subvolume ranks remain fixed while both rank-OLS projections are refitted in each centre-coordinate spatial-block bootstrap replicate (2,000 accepted replicates per target).

The specified main associations are lambda5 with tail imbalance, and lambda6 with outer-tail probability and Q95-Q05. `conditional_rank_association_diagnostics.csv`, `conditional_residual_balance.csv`, and `numerical_support.csv` retain the adjustment and support diagnostics. `conditional_ecdf.csv`, `conditional_ecdf_contrasts_bootstrap.csv`, `conditional_tail_summary.csv`, and `conditional_tail_contrasts_bootstrap.csv` supply ECDF and high-Q80-minus-low-Q20 summaries; `illustrative_conditional_pairs.csv` records illustrative windows.

These are conditional associations among descriptors computed from the same local marginal distributions. The unexplained rank variance is a descriptive non-redundancy measure. ECDF groups and selected pairs are descriptive illustrations; they do not establish external predictive performance or independent spatial-connectivity information. All five original supplement PNGs are included.

## Figures and compact plotting data

The six direct-statistics PNGs are in `sections_4_4_to_4_6/figures/`; the three directional-template PNGs are in the Section 4.4 component's `figures/`. Original descriptive filenames are retained. `report_artifact_contract.csv` records the direct figure/table filters, source hashes, and color scales.

`figure_inputs/direct_representative_z50.csv` contains 27,264 rows for z45_56 at z=50, primary 5x5x3 window: Section 4.5 uniform plus the three Section 4.6 weightings. It supplies the representative maps and Section 4.5 Spearman correlation matrix. Retain the section distinction even where uniform appears in both sections, and retain `analysis_valid` and the recorded support scope.

`figure_inputs/template_2d_original_and_rotated.csv` contains all 3,888 original/rotated rows across the recorded datasets, properties, templates, and evaluators. The surrogate-summary CSV supplies the 2D envelopes; the complete 3D response table is included in the component. These selections support the plots but do not contain the complete bootstrap resampling population or surrogate ensemble.

## Integrity and reproduction

- `export_selection.json` lists the exact-copy selection. `copied_files_manifest.csv` gives source-relative and published paths, byte counts, and SHA-256 for all 234 copied files. All 174 copied direct-run files and 18 supplement outputs also matched their source archive manifests during curation.
- `figure_inputs/derived_manifest.json` records selections, source hashes, source/selected row counts, and output hashes. Both large source CSV hashes were checked against the source archive before extraction.
- Section 4.3 `sha256_manifest.csv` is the original execution-provenance manifest for code, build inputs, binaries, and inputs. The supplement manifest also lists large omitted outputs. Neither is an inventory of this compact package. The curated copy manifest covers the files distributed here.
- Source JSON, INI, report, and manifest files retain recorded machine paths. They identify original inputs and outputs; adapt input/build/output paths in the supplied repository runners for a new run. A copied configuration is provenance, not a portable execution command.

From the repository root, verify the distributed bytes with:

```text
python Paper-data/tools/verify_current_results.py
```

Recreate the compact CSVs from the complete source archive into a new output directory with:

```text
python Paper-data/tools/export_figure_inputs.py --source-run <complete-direct-run-directory> --output-dir <new-output-directory>
```

Full regeneration uses Stanford VI-E porosity/facies files from the [SCRF repository](https://github.com/SCRFpublic/Stanford-VI-E), `standalone/experiments/stanford_vie_preprocessing/`, the Section 4.3 benchmark, `standalone/experiments/run_direct_sections_44_46.py`, and `standalone/experiments/r1_5_incremental_value_20260924/conditional_high_order_window_analysis.py`. Use the matching frozen protocols. The preprocessing manifest records source/output hashes, axes, units, and fixtures, including the `poro.inc` coarse-grid check. Current slabs are z=15-26, 45-56, and 75-86; the benchmark uses the separately recorded z=62 slice. Recorded source paths use `WithoutSGeMS`; the public implementation is under `standalone/`.

Omitted full data include per-cell maps/support masks, the 510,097,911-byte direct analysis table, the 89,790,107-byte 2D template response table, full H2 acceptance/diagnostic histories, generated surrogate inputs, logs, binaries, the 112,089,799-byte conditional residual table, and the 12,433,352-byte conditional association bootstrap table. Recomputing all intervals and surrogate tests requires regenerating those full inputs. Recorded PASS states describe the archived execution and do not certify a new build or timing run.

Stanford VI-E retains its original attribution and repository terms; this curation does not assign new rights to source data.
