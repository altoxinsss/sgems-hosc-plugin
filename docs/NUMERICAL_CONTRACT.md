# HOSC SGeMS Numerical Contract

## Supported inputs

`HOSC_WINDOW` and the formal `HOSC_TEMPLATE` paths accept an informed `Cartesian_grid` only. Template and window coordinates are derived from the grid node index and its dimensions. Ordinary properties named `x`, `y`, or `z` are never interpreted as geometry. A zero property value remains an observation; uninformed and non-finite values are excluded.

The moving-window route requires every member of a window to be informed and finite. An incomplete window has no primary result. Its diagnostics record the available support and a numeric-valid value of zero.

## Window statistics

`Target_Button=Cumulant` supports orders 1 through 8. For order 1 the result is the requested first cumulant. When `center_cumulant_input_for_order_ge_2=true` (the default), cumulant inputs of order 2 and above are translated by their window mean before Bell, AVL, or the shared weighted recurrence is evaluated. This control does not change the `Moment` route. With soft weighting, weights are computed from the original values, normalized, and then held fixed while samples are translated.

The primary result property stores raw kappa_r. When `output_standardized_cumulant=true`, the plugin also writes `<output>_kappa2`, `<output>_lambda`, and `<output>_lambda_valid`, where `lambda_r = kappa_r / kappa_2^(r/2)`. The numerator and denominator use the same complete support and normalized weights. `lambda_valid` is zero when kappa_2 is non-finite or no greater than `variance_epsilon`; the lambda property is then uninformed. Raw kappa_r is retained whenever it is numerically valid.

Primary-statistic inversion and legacy first-moment centering are fixed to false in the backend and have no UI controls. Old parameter files may omit `invert_window_output` and `zero_first_order_moment`, or specify false; explicit true values are rejected before calculation. Stable cumulant translation is controlled only by `center_cumulant_input_for_order_ge_2`.

## Weighting

The supported moving-window modes are `uniform`, `coordinate_gaussian`, and `soft_high_porosity`. Non-uniform modes require `enable_weighted_moment=true`. Weight normalization is fixed to true in the backend and has no UI control. Old parameter files may omit `normalize_weights` or specify true; an explicit false value is rejected before calculation, including when weighting is disabled.

`coordinate_gaussian` always uses the window center. For every window, its geometric center and per-axis half-window extent define the relative coordinate Gaussian. Singleton axes contribute zero normalized distance. The fixed center policy has no UI control. Old parameter files may omit `weight_center_mode` or specify `window_center`; conflicting explicit values are rejected before calculation.

`soft_high_porosity` uses `(1-eta) + eta * sigmoid((value-threshold)/transition_width)` with the configured threshold, transition width, and strength.

## Template sampling

The template input boundary is fixed at zero-based Cartesian node indices; missing source nodes do not alter that basis. After direction reordering and before calculation begins, the plugin freezes the displayed direction count, result-row columns, lag steps, and output dimensions. Result properties use SGeMS x-fastest node order, with every result row placed from this frozen layout rather than from the internal lag-enumeration order or calculator-owned `Keys` arrays. Three selected display directions remain a three-direction layout even when the third output dimension is one.

All public `HOSC_TEMPLATE` calculation routes are unweighted, and the interface has no weighting tab. Old parameter files may omit `Weighted_Enable` or specify `Disable`; enabling it is rejected before either legacy or exact sampling. This also applies when a legacy file selects a uniform weight mode.

`HOSC_TEMPLATE` retains `template_sampling=legacy_search` only for 2D Cartesian grids. Three-dimensional template calculations require `template_sampling=exhaustive_integer_offsets`. The exact route is a cumulant path: it requires one signed nonzero integer vector per configured direction in `integer_direction_vectors`, constructs offsets from the selected lag cell, and evaluates all complete translations of that fixed template. `LAG` is a positive cell step, individual selected lag indices may be zero, and the total joint cumulant order (origin plus all direction occurrences) is at most eight. A lag with no complete translations is uninformed rather than zero.

## Output and audit rules

Numeric validity, lambda validity, and float display validity are separate. Before a double result is written to an SGeMS property, it must be finite and within the finite `float` range. Float-overflow values remain available only in a requested audit sidecar and are uninformed in the host property.

`write_audit_sidecar=true` requires an existing `audit_output_directory`. Each run creates a unique `HOSC_WINDOW_*` directory containing `resolved_config.txt` and `window_results.csv`. Audit schema `hosc_window_audit.v3` records the build identity, source grid/property, all active weight parameters, centering and output switches, and resolved geometry. The CSV records double raw values, kappa2, lambda, support, ESS, all three validity states, whether lambda was requested, window start, window center, and a reason. `lambda_not_requested` is distinct from a low-variance lambda failure. This sidecar is the comparison input for host numerical acceptance; it is not a replacement for the saved SGeMS project.
