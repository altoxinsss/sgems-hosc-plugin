param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = "Stop"
$pluginPath = Join-Path $RepoRoot "MoveWindow\HOSCWindow\Plugin.cpp"
$corePath = Join-Path $RepoRoot "shared\WeightedMomentCore.cpp"
$uiPath = Join-Path $RepoRoot "ui\HOSC_WINDOW.ui"
$plugin = Get-Content -LiteralPath $pluginPath -Raw
$core = Get-Content -LiteralPath $corePath -Raw
$ui = Get-Content -LiteralPath $uiPath -Raw

function Assert-Match([string]$Text, [string]$Pattern, [string]$Message) {
    if ($Text -notmatch $Pattern) { throw $Message }
}

Assert-Match $plugin 'center_cumulant_input_for_order_ge_2' "The stable cumulant-centering control is missing."
Assert-Match $plugin 'output_standardized_cumulant' "The standardized cumulant output control is missing."
Assert-Match $plugin 'variance_epsilon' "The standardized cumulant variance threshold is missing."
Assert-Match $plugin 'G_CalculateTarget\s*!=\s*4' "Standardized output is not restricted to Cumulant."
Assert-Match $plugin 'Parameter center mode is not supported' "Coordinate Gaussian parameter-center mode is not explicitly rejected."
Assert-Match $plugin 'source_nx_\s*==\s*1\s*\|\|\s*source_ny_\s*==\s*1' "The unsupported XZ/YZ-grid guard is missing."
Assert-Match $plugin 'source_present_\[\(size_t\)\(source_index\s*-\s*1\)\]\s*==\s*0' "Missing source cells are not excluded from complete windows."
Assert-Match $plugin 'is_float_representable\(value\)' "Float32 display validity is not checked before host output."
Assert-Match $plugin 'write_window_audit_sidecar' "Double-precision audit sidecar writing is missing."
Assert-Match $plugin 'hosc_window_audit\.v3' "The audit sidecar schema is not v3."
Assert-Match $plugin 'hosc_window_run_status\.v1' "The per-run audit status record is missing."
Assert-Match $plugin 'enable_percentile_clipping' "The audit sidecar does not record the resolved clipping switch."
Assert-Match $plugin 'percentile_vmin' "The audit sidecar does not record the resolved lower clipping percentile."
Assert-Match $plugin 'percentile_vmax' "The audit sidecar does not record the resolved upper clipping percentile."
Assert-Match $plugin 'window_params\.call\.returned' "The parameter-function return boundary trace is missing."
Assert-Match $plugin 'window_ready\.host_log\.begin' "The host window-ready log boundary trace is missing."
Assert-Match $plugin 'window_params\.body\.end\.locals_destroy_next' "The parameter local-destruction boundary trace is missing."
Assert-Match $plugin 'lambda_not_requested' "The audit sidecar cannot distinguish disabled lambda output from a variance failure."
Assert-Match $plugin 'Inverting the primary statistic is not supported' "The primary-statistic inversion guard is missing."
Assert-Match $plugin 'Legacy zero-first-order centering is not supported' "The unsupported legacy centering guard is missing."
Assert-Match $plugin 'is_float_representable\(clipped_values' "The clipped-property float32 range guard is missing."
Assert-Match $plugin 'source_grid_name_' "The audit sidecar does not receive source-grid identity."
Assert-Match $plugin 'source_property_name_' "The audit sidecar does not receive source-property identity."
if ($plugin -match 'tempData\.txt') { throw "The legacy fixed tempData.txt export remains in the active window plugin." }
Assert-Match $core 'configure_relative_window_coordinate_gaussian' "The compiled shared core lacks relative-window coordinate Gaussian support."
Assert-Match $core 'weighted_univariate_cumulant' "The compiled shared core lacks the univariate cumulant route."
Assert-Match $ui '<cstring>center_cumulant_input_for_order_ge_2</cstring>' "The UI has no cumulant-centering control."
Assert-Match $ui '<cstring>output_standardized_cumulant</cstring>' "The UI has no standardized-cumulant control."
Assert-Match $ui '<cstring>variance_epsilon</cstring>' "The UI has no variance-epsilon control."
Assert-Match $ui '<cstring>write_audit_sidecar</cstring>' "The UI has no audit-sidecar control."
Assert-Match $ui 'Invert primary statistic \(unsupported\)' "The UI still presents primary-statistic inversion as supported."
Assert-Match $ui 'Legacy centering \(unsupported\)' "The UI still presents legacy centering as supported."
Assert-Match $plugin 'hosc_validation::FeedbackGuard feedback\(error_mesgs, "hosc_window_validation_summary"\)' "The moving-window plugin does not publish its input errors to the form."
Assert-Match $ui '<cstring>hosc_window_validation_summary</cstring>' "The moving-window UI has no field-specific validation summary."

$source = @(75, 100, 12)
$window = @(5, 5, 3)
$step = @(1, 1, 1)
$output = for ($i = 0; $i -lt 3; $i++) { 1 + [math]::Floor(($source[$i] - $window[$i]) / $step[$i]) }
if (($output -join ',') -ne '71,96,10') { throw "Full-window dimensions are incorrect: $($output -join 'x')" }

$center = for ($i = 0; $i -lt 3; $i++) { 0.5 * ($window[$i] - 1) }
if (($center -join ',') -ne '2,2,1') { throw "Window-center origin is incorrect: $($center -join ',')" }

Write-Host "Moving-window static contract verification passed against current source and UI. This does not execute the host or numerical core."
