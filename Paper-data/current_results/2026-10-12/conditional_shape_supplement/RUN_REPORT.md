# Conditional High-Order Window Analysis

## Scope

This post-processing analysis uses only the frozen run-local porosity grids and uniform 5x5x3 Section 4.5 numerical maps. Channel labels, facies validity, prediction models, and C++ reruns are not used.

For `lambda5`, the rank projection conditions on mean, log variance, high-porosity proportion, lambda3, and lambda4. For `lambda6`, it additionally conditions on lambda5. The Q20 and Q80 groups are defined from the resulting rank-OLS residuals within each subvolume.

## Numerical Support

| dataset_id | structural_rows | numerical_windows | window_samples | reference_variance | variance_epsilon | minimum_variance_retained | minimum_ess | maximum_ess | formal_analysis_variance_flag_disagreement |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z15_26 | 68160 | 68160 | 75 | 0.007518 | 7.518e-11 | 5.965e-06 | 75 | 75 | 0 |
| z45_56 | 68160 | 68160 | 75 | 0.006645 | 6.645e-11 | 8.977e-06 | 75 | 75 | 0 |
| z75_86 | 68160 | 68160 | 75 | 0.006065 | 6.065e-11 | 1.718e-05 | 75 | 75 | 0 |

Retained windows are required independently to have finite numerical fields, `cpp_valid=true`, effective sample size 75, and `kappa2 > variance_epsilon`. The accompanying `numerical_support.csv` records the exclusion counts for each gate.

## Conditional Rank Variation

| dataset_id | target | controls | r_squared | unexplained_rank_variance | low_q20_count | high_q80_count |
| --- | --- | --- | --- | --- | --- | --- |
| z15_26 | lambda5 | mean;log_variance;high_fraction;lambda3;lambda4 | 0.1105 | 0.8895 | 13632 | 13632 |
| z45_56 | lambda5 | mean;log_variance;high_fraction;lambda3;lambda4 | 0.2267 | 0.7733 | 13632 | 13632 |
| z75_86 | lambda5 | mean;log_variance;high_fraction;lambda3;lambda4 | 0.2729 | 0.7271 | 13632 | 13632 |
| z15_26 | lambda6 | mean;log_variance;high_fraction;lambda3;lambda4;lambda5 | 0.2451 | 0.7549 | 13632 | 13632 |
| z45_56 | lambda6 | mean;log_variance;high_fraction;lambda3;lambda4;lambda5 | 0.4117 | 0.5883 | 13632 | 13632 |
| z75_86 | lambda6 | mean;log_variance;high_fraction;lambda3;lambda4;lambda5 | 0.6697 | 0.3303 | 13632 | 13632 |

The unexplained rank variance is a descriptive non-redundancy measure. It is not a prediction score or an external geological-validation metric.

## Fixed-Rank Conditional Associations

For each empirical-distribution descriptor, its complete-subvolume rank and the cumulant rank are separately projected on the same lower-order rank controls. The reported association is the correlation of those two residuals. Within every centre-coordinate spatial-block bootstrap replicate, both projections are re-fit from the resampled block sufficient statistics; the full-subvolume rank transformation remains fixed.

| target | metric | conditional_rank_association | ci95_low | ci95_high | conditional_residual_slope | bootstrap_replicates |
| --- | --- | --- | --- | --- | --- | --- |
| lambda5 | tail_imbalance | -0.649 | -0.6577 | -0.638 | -0.3497 | 2000 |
| lambda6 | outer_tail_probability | -0.4757 | -0.4862 | -0.4607 | -0.5895 | 2000 |
| lambda6 | outer_span_q05_q95 | -0.5788 | -0.5969 | -0.5555 | -0.616 | 2000 |

The three rows above were specified before this analysis: lambda5 with tail imbalance, and lambda6 with outer-tail probability and Q95-Q05. `conditional_rank_associations.csv` contains the complete descriptor matrix. The local empirical-distribution descriptors and cumulants are calculated from the same windows, so this is a conditional shape-association analysis rather than an external prediction experiment.

| dataset_id | target | metric | target_rank_r_squared | metric_rank_r_squared | target_max_abs_residual_control_correlation | metric_max_abs_residual_control_correlation | centre_coordinate_block_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| z15_26 | lambda5 | tail_imbalance | 0.1105 | 0.7623 | 7.494e-16 | 2.151e-15 | 144 |
| z15_26 | lambda6 | outer_tail_probability | 0.2451 | 0.1381 | 4.25e-16 | 4.136e-16 | 144 |
| z15_26 | lambda6 | outer_span_q05_q95 | 0.2451 | 0.3434 | 4.25e-16 | 8.197e-16 | 144 |
| z45_56 | lambda5 | tail_imbalance | 0.2267 | 0.7829 | 3.836e-16 | 1.042e-15 | 144 |
| z45_56 | lambda6 | outer_tail_probability | 0.4117 | 0.1103 | 3.996e-16 | 6.005e-16 | 144 |
| z45_56 | lambda6 | outer_span_q05_q95 | 0.4117 | 0.4252 | 3.996e-16 | 8.901e-16 | 144 |
| z75_86 | lambda5 | tail_imbalance | 0.2729 | 0.7677 | 2.954e-15 | 3.94e-15 | 144 |
| z75_86 | lambda6 | outer_tail_probability | 0.6697 | 0.2631 | 1.589e-15 | 1.341e-15 | 144 |
| z75_86 | lambda6 | outer_span_q05_q95 | 0.6697 | 0.4207 | 1.589e-15 | 7.174e-16 | 144 |

## Block-Bootstrap Local-Distribution Contrasts

| target | metric | contrast_high_q80_minus_low_q20 | ci95_low | ci95_high | bootstrap_replicates |
| --- | --- | --- | --- | --- | --- |
| lambda5 | tail_imbalance | -0.1246 | -0.1301 | -0.1191 | 2000 |
| lambda6 | outer_tail_probability | -0.08388 | -0.08631 | -0.08143 | 2000 |
| lambda6 | outer_span_q05_q95 | -0.5793 | -0.5993 | -0.5595 | 2000 |

Each contrast is an equal-subvolume average. Spatial blocks, rather than individual samples, were resampled by window centre coordinate. The ECDF figures and Q20/Q80 contrasts remain descriptive group summaries; the fixed-rank conditional associations above provide the lower-order-adjusted quantitative comparison. Outer-tail summaries are compact descriptors rather than a claim that a single cumulant represents only one tail property.

## Illustrative Window Pairs

| dataset_id | target | selected | nearest_candidates_examined_per_high_anchor | high_anchor_coverage_within_examined_neighbours | eligible_nonoverlapping_pair_count | high_x | high_y | high_z | low_x | low_y | low_z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z15_26 | lambda5 | true | 8 | 0.9939 | 104469 | 38 | 57 | 20 | 36 | 63 | 18 |
| z45_56 | lambda5 | true | 8 | 0.9994 | 106033 | 73 | 63 | 55 | 38 | 97 | 53 |
| z75_86 | lambda5 | true | 8 | 0.8425 | 81484 | 13 | 47 | 78 | 63 | 54 | 82 |
| z15_26 | lambda6 | true | 8 | 0.9999 | 105150 | 69 | 17 | 17 | 36 | 88 | 21 |
| z45_56 | lambda6 | true | 8 | 0.8921 | 92961 | 42 | 24 | 49 | 57 | 67 | 52 |
| z75_86 | lambda6 | true | 8 | 0.9905 | 105158 | 54 | 97 | 82 | 19 | 54 | 83 |

The selected pairs are traceable illustrations. The quantitative evidence is based on all residual-group windows and the spatial-block bootstrap, not on a selected pair.

## Bootstrap Diagnostics

```json
{
  "fixed_rank_conditional_associations": {
    "lambda5": {
      "accepted": 2000,
      "attempts": 2000
    },
    "lambda6": {
      "accepted": 2000,
      "attempts": 2000
    }
  },
  "residual_group_contrasts": {
    "lambda5": {
      "accepted": 2000,
      "attempts": 2000
    },
    "lambda6": {
      "accepted": 2000,
      "attempts": 2000
    }
  }
}
```

## Source Hashes

```json
{
  "analysis_utc": "2026-10-11T23:11:57.565903+00:00",
  "conditional_association_protocol": {
    "bootstrap_refits_rank_ols": true,
    "bootstrap_resampling_unit": "centre_coordinate_spatial_block",
    "primary_metrics": {
      "lambda5": [
        "tail_imbalance"
      ],
      "lambda6": [
        "outer_tail_probability",
        "outer_span_q05_q95"
      ]
    },
    "rank_transform_scope": "complete_subvolume_fixed"
  },
  "direct_analysis_values": {
    "run_relative_path": "direct_analysis_values.csv",
    "sha256": "25af0f82e9e8eb5ab45238ac868e47ffe10a050fae5425aa209d12a45629b91e"
  },
  "formal_evidence_manifest": {
    "run_relative_path": "evidence_manifest.json",
    "sha256": "01bc348b461c8dbfd47d1b45762167a6071b5363df0000a3a4508137ec09b678"
  },
  "formal_gate": {
    "run_relative_path": "gate.json",
    "sha256": "bbddf101e2e99e6ffda586d5ea590f059db5772fdd359f37183296ba23a3ef8c"
  },
  "formal_run_directory": "D:\\Development\\DLL\\DLL\\WithoutSGeMS\\output\\direct_4_4_4_6_20261011\\20261012_031025\\run_20261011_191032",
  "formal_sha256_inventory": {
    "run_relative_path": "sha256_manifest.csv",
    "sha256": "4be44c0d125a76fe97a3d27a825112b400d656c571a9c53bac222bc5fb7b0013"
  },
  "frozen_common_protocol": {
    "run_relative_path": "configs/direct_common_protocol.json",
    "sha256": "29a28bf63507cea7ba7dd0b054c3fbe3b38ac303df6f623a938de16beea9e2fe"
  },
  "frozen_section45_protocol": {
    "run_relative_path": "configs/direct_section45.json",
    "sha256": "8ffafe1e27e82d42ac4bc2d9069ae3061cacfb4bdf853fc024f6f0311350fb84"
  },
  "high_porosity_threshold": 0.25,
  "preprocessing_manifest": {
    "run_relative_path": "preprocessing/preprocessing_manifest.json",
    "sha256": "83f6945a893c6870070e5b8df0a646fed20c8484113a9955819e19c197a3c439"
  },
  "required_effective_sample_size": 75,
  "run_local_porosity_grids": {
    "z15_26": {
      "reference_variance": 0.007517608758404953,
      "run_relative_path": "preprocessing/z15_26/stanford_vie_z15_26_porosity.grid",
      "sha256": "209f78e6edd798de3bc8079ffcc7b2df0d10eb15d3407b9622ca191d31217142",
      "variance_epsilon": 7.517608758404952e-11
    },
    "z45_56": {
      "reference_variance": 0.006645234264241945,
      "run_relative_path": "preprocessing/z45_56/stanford_vie_z45_56_porosity.grid",
      "sha256": "227fa518d2d743a7ab13c4a3180b9bbc915229de958737db009f323d5d182974",
      "variance_epsilon": 6.645234264241945e-11
    },
    "z75_86": {
      "reference_variance": 0.006065063612806557,
      "run_relative_path": "preprocessing/z75_86/stanford_vie_z75_86_porosity.grid",
      "sha256": "8e2e50ad2d4d126e0fe48998711105b7c4743f54fade908e6575d424564bb1d4",
      "variance_epsilon": 6.065063612806557e-11
    }
  },
  "run_manifest": {
    "run_relative_path": "run_manifest.json",
    "sha256": "da0db4358af6b5ec773e9a09fba56401fbf4572a58de7372748884ca1a1cbb5c"
  },
  "script": "D:\\Development\\DLL\\DLL\\WithoutSGeMS\\experiments\\r1_5_incremental_value_20260924\\conditional_high_order_window_analysis.py",
  "script_sha256": "043d066f163c4220ed1e12a16729caba0becb0c76c832c476f454038f0d44bc2",
  "spatial_block_shape_xyz": [
    8,
    13,
    6
  ],
  "uniform_window": [
    5,
    5,
    3
  ]
}
```
