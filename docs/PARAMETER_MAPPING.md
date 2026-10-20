# HOSC SGeMS Parameter Mapping

## HOSC_WINDOW

`HOSC_WINDOW` operates only on an informed `Cartesian_grid`. It reads the
selected property's grid-node order directly; ordinary properties named `x`,
`y`, or `z` do not define geometry.

| UI/XML key | Meaning | Contract |
|---|---|---|
| `Source_Grid.value` | source grid | Required Cartesian grid |
| `Source_Property.value` | source property | Required and read with the informed mask |
| `props.value` | output property name | Required |
| `Target_Button.value` | calculation target | `Moment` or `Cumulant` |
| `Algorithm_Button.value` | unweighted evaluator | `Bell` or `AVL` |
| `windowl/windoww/windowh.value` | window node counts | Positive integers; complete window required |
| `stepx/stepy/stepz.value` | window movement | Positive in every non-singleton axis; `stepz=0` only for XY single-layer input |
| `moment_order.value` | requested order | Integer 1 through 8 |
| `output_diagnostics.value` | output diagnostics | Boolean; defaults to true |
| `center_cumulant_input_for_order_ge_2.value` | stable cumulant translation | Boolean; default true; applies only to Cumulant order >=2 |
| `output_standardized_cumulant.value` | derive lambda output | Boolean; Cumulant only |
| `variance_epsilon.value` | minimum kappa2 for lambda | Finite value >=0; default `1e-12` |
| `invert_window_output.value` | invert primary output | Only `false` is supported; `true` is rejected to preserve the kappa/lambda definition |
| `zero_first_order_moment.value` | legacy weighted-moment centering | Only `false` is supported; `true` is rejected; use `center_cumulant_input_for_order_ge_2` for cumulants |
| `write_audit_sidecar.value` | write unique double audit sidecar | Boolean; default false |
| `audit_output_directory.value` | audit base directory | Required when `write_audit_sidecar=true`; must be an existing writable directory |
| `enable_weighted_moment.value` | weight route switch | Boolean |
| `weight_mode.value` | weighting scheme | `uniform`, `coordinate_gaussian`, or `soft_high_porosity` |
| `normalize_weights.value` | normalized weighted statistic | Must be true when weighting is enabled |
| `weight_bandwidth.value` | coordinate Gaussian sigma | Finite and greater than zero |
| `weight_center_mode.value` | coordinate center choice | `window_center` by default; `parameter` is rejected |
| `high_porosity_threshold.value` | soft-weight threshold | Finite value; UI default 0.25 |
| `soft_transition_width.value` | soft-weight transition width | Finite and greater than zero; UI default 0.02 |
| `soft_weight_strength.value` | soft-weight strength | Finite value in [0, 1]; UI default 0.8 |

For `coordinate_gaussian`, the implementation resolves a separate geometric
center for each output window. Per-axis scale is half the window extent in
that axis. Singleton axes contribute zero normalized distance. Parameter-center
mode is not supported in this build.

The weighted route uses the shared weighted recurrence. The `Bell`/`AVL`
selection describes the unweighted cumulant evaluator and is not represented
as two independent weighted implementations.

When diagnostics are enabled, the output grid also receives properties named
`<output>_support`, `<output>_ess`, `<output>_numeric_valid`, and
`<output>_display_valid`. The primary output remains uninformed for an
incomplete window, non-positive/non-finite weight support, a non-finite final
statistic, or float32 display overflow. Standardized output also adds
`<output>_kappa2`, `<output>_lambda`, and `<output>_lambda_valid`.

## HOSC_TEMPLATE

| UI/XML key | Meaning | Contract |
|---|---|---|
| `template_sampling.value` | template construction route | `legacy_search` for 2D only, or `exhaustive_integer_offsets` for 2D/3D |
| `integer_direction_vectors.value` | signed integer offsets | Exact mode: one nonzero `dx,dy[,dz]` vector per configured direction, separated by semicolons |
| `NOD/NLAG/LAG/LTOL/...` | direction records | All fields, including LTOL, move together when directions are rearranged for output; `LAG>0`; origin plus occurrence counts must be at most order 8 |

The published `HOSC_TEMPLATE` UI exposes uniform weighting only. Its legacy
weighted implementation remains source-compatible with historic XML files but
is not a formal UI route. `exhaustive_integer_offsets` accepts unweighted
Cumulant only, applies complete-support sampling, uses fixed zero-based source
coordinates, and leaves a lag without complete translations uninformed rather
than exporting zero. Result lag rows are written to the x-fastest SGeMS output
node determined by their displayed lag indices. The selected display-direction
count is retained even when a displayed axis has `NLAG=1`; hidden directions
are represented only by the fixed cross-section coordinates. Template and
moving-window computations remain distinct: the template algorithm evaluates
a fixed multi-location structure, whereas `HOSC_WINDOW` calculates local
single-variable statistics at moving window centers.
