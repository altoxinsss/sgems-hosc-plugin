# SGeMS parameter correspondence

The standalone programs and SGeMS plugins provide template joint cumulants and
moving-window marginal cumulants. This document maps the current supported
settings; it does not certify SGeMS host runtime stability. Plugin interfaces
are specified in [PARAMETER_MAPPING.md](../docs/PARAMETER_MAPPING.md), with
host checks in [MANUAL_TESTING.md](../docs/MANUAL_TESTING.md).

## Moving-window settings

Set the following explicitly in the CLI INI. These are a common comparison
profile, not a statement that every CLI default matches the plugin default.

```ini
[general]
calculate_target=Cumulant
calculate_algorithm=Bell
; Set these to the source geometry used in the comparison.
source_origin=0,0,0
source_spacing=1,1,1
; 0 requests automatic detection; set 1, 2, or 3 when needed.
coordinate_dimensions=0

[sliding_window]
moment_cumulant_order=4
window_size=3,3,3
move_step=1,1,1
require_full_window=true
require_complete_window=true
require_anchor_present=true
skip_missing_samples=false
output_at_window_center=true
center_cumulant_input_for_order_ge_2=true
invert_window_output=false
enable_percentile_clipping=false

[weighted_moment]
enable_weighted_moment=false
weight_mode=uniform
enable_coordinate_gaussian=false
enable_soft_high_porosity=false
normalize_weights=true
zero_first_order_moment=false
weight_center_mode=window_center
bandwidth=1.0
min_weight=0
high_porosity_threshold=0.25
soft_transition_width=0.02
soft_weight_strength=0.8

[runtime]
write_window_support_sidecar=true
```

Supply `source_grid`, `source_property`, `target_property`, and `output_dir`
for the actual input and destination. Paths are resolved relative to the INI.
Use the [current weighting fragments](configs/weighting/) to change schemes.
`bandwidth` is the CLI equivalent of the plugin's `weight_bandwidth` key.
The current SGeMS window interface requires normalized weights and rejects
`zero_first_order_moment=true`, output inversion, and a parameter-fixed
coordinate centre. Set `window_center` explicitly in the CLI: its retained
default is `parameter`, whereas the current plugin uses `window_center`.
The CLI centering default is also `false`, so set the centering flag explicitly.

Only complete in-bounds windows with complete finite data support contribute
to this profile. Missing values are represented by the grid's informed state,
not by the numeric value zero. The plugin leaves invalid output nodes
uninformed; the CLI omits their text rows and supplies a window-support CSV.
These file/property representations must be matched by coordinates and masks
before comparing values.

## Geometry and outputs

The source must describe a regular Cartesian lattice. Columnar CLI input may
contain one to three coordinate columns followed by named value columns.
`source_property` selects a value-column header case-insensitively; for
`x y porosity permeability`, set `coordinate_dimensions=2` when selecting
`permeability`. The SGeMS plugin instead obtains node geometry from its
`Cartesian_grid`; ordinary properties named `x`, `y`, or `z` are not coordinates.

The CLI metadata retains the floating-point centre of the first window as
`target_origin`. Text-grid row coordinates are integer lattice indices.
Source origin, spacing, output steps, and index conventions must therefore be
checked when displaying or comparing maps.

The plugin can write raw cumulants, standardized `lambda`, validity properties,
and a double-precision audit sidecar. The standalone direct-statistics workflow
derives standardized cumulants and its common-support diagnostics from the
CLI results. The two interfaces do not produce identical file sets. Compare
raw statistics on the same support and normalization; plugin display properties
are float32, whereas its audit sidecar retains double values.

## Directional templates

The current formal template comparison uses unweighted Cumulant mode with
`template_sampling=exhaustive_integer_offsets`. Supply the same signed integer
direction vectors, lag indices, occurrence counts, source geometry, and
complete-support rules in each interface. The SGeMS template UI exposes uniform
weighting only. Template lag-response maps and moving-window maps have different
coordinate meanings and must not be compared as interchangeable outputs.

For build commands and executable locations, use the
[standalone README](README.md). For complete scientific runs, use the
[reproduction guide](experiments/README.md).
