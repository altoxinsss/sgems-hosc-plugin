# Current weighting schemes

Formula version: `2026-10-10.v2`.

The formal moving-window workflow uses `uniform`, `coordinate_gaussian`, and
`soft_high_porosity`. Section 4.5 uses uniform weights; Section 4.6 compares
the two non-uniform schemes with the same uniform baseline on common support.
The formal directional-template workflow uses uniform weights.

## Normalization and cumulants

For one complete window, let `g_i` be the nonnegative raw weight and define
`p_i = g_i / sum_j(g_j)`. Moments and cumulants use these normalized empirical
weights. The effective sample size is `ESS = 1 / sum_i(p_i^2)`.

Use `normalize_weights=true`, `min_weight=0`, and
`zero_first_order_moment=false`. For cumulants of order two and above, set
`center_cumulant_input_for_order_ge_2=true` in `[sliding_window]`.
This numerical translation preserves cumulants; value weights are calculated
from the original observations and held fixed during translation.

## Three schemes

| Mode | Raw weight | Formal primary settings |
|---|---|---|
| `uniform` | `g_i = 1` | No value or location preference |
| `coordinate_gaussian` | `g_i = exp[-sum_j(u_ij^2) / (2 sigma_s^2)]` | `bandwidth=1.0`, `weight_center_mode=window_center` |
| `soft_high_porosity` | `g_i = (1 - eta_v) + eta_v sigmoid((z_i - t) / b)` | `t=0.25`, `b=0.02`, `eta_v=0.8` |

For coordinate weights, `u_ij` is the offset from the geometric window centre
divided by that axis's half-window extent. Singleton axes contribute zero.
Using physical coordinates and the corresponding physical half-extent cancels
the grid spacing, so `sigma_s` is dimensionless. These weights emphasize the
window centre and do not encode a directional multi-point template.

For value weights, `sigmoid(a) = 1 / (1 + exp(-a))`.
The INI keys are `high_porosity_threshold=t`, `soft_transition_width=b`, and
`soft_weight_strength=eta_v`. The weight increases smoothly with porosity;
with `eta_v=0.8`, its limiting raw values are 0.2 and 1. This scheme emphasizes
a value range and is not a sampling-density correction.

Primary and sensitivity settings are frozen in
[direct_section46.json](experiments/configs/direct_section46.json).
The three current [INI fragments](configs/weighting/) supply matching examples.

## Selecting the route

`weight_mode` selects a scheme. For a weighted CLI evaluation, also set
`enable_weighted_moment=true`; the unweighted Bell/AVL baseline uses `false`.
If explicit base Boolean switches are supplied, they override `weight_mode`:
set `enable_coordinate_gaussian` or `enable_soft_high_porosity` to `true` for
the corresponding mode, and the other to `false`. Both `false` select uniform.
More than one enabled base switch is an error.

The shared standalone core retains compatibility branches for older inputs,
but the formal experiment generator rejects fixed value-Gaussian, adaptive
upper-tail, and joint-geometric weighting. They are not current experiment
schemes or published example configurations.

## Reproducibility

The CLI metadata records `weighting_formula_version`, `weighting_scheme_id`,
normalization, centering, and resolved weighting parameters. Keep the generated
INI, output metadata, executable and source hashes, dataset hashes, and run
gate together. The Python reference is
`experiments/r1_5_incremental_value_20260924/weighting.py`; the formal CLI
configuration writer is `cli_runner.py` in the same directory.

See [SGeMS parameter correspondence](SGEMS_COMPATIBILITY.md) for explicit
settings and the differences between CLI files and plugin properties.
