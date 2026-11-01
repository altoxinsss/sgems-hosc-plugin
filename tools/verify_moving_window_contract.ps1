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
Assert-Match $plugin 'Coordinate Gaussian weights always use the window center' "The fixed coordinate-Gaussian center policy is missing."
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
Assert-Match $plugin 'The primary statistic is never inverted' "The primary-statistic inversion guard is missing."
Assert-Match $plugin 'Legacy zero-first-order centering is disabled' "The legacy centering guard is missing."
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
foreach ($fixedField in @('normalize_weights', 'zero_first_order_moment', 'invert_window_output', 'weight_center_mode')) {
    if ($ui -match ('<cstring>' + [regex]::Escape($fixedField) + '</cstring>')) {
        throw "The fixed backend setting $fixedField is still exposed in the UI."
    }
}
if ($ui -match 'Normalize weights|Legacy centering|Invert primary statistic|Coordinate center mode') {
    throw "A removed fixed-setting label remains in the moving-window UI."
}

# Legacy PAR values are validated separately from the fixed execution state.
# In particular, normalization=false must not be accepted merely because the
# current run has weighting disabled.
Assert-Match $plugin 'weightedOptions\.normalize_weights\s*=\s*true;' "Window normalization is not fixed to true."
Assert-Match $plugin 'weightedOptions\.zero_first_order_moment\s*=\s*false;' "Legacy first-moment centering is not fixed to false."
Assert-Match $plugin 'invertWindowOutput\s*=\s*false;' "Primary-statistic inversion is not fixed to false."
Assert-Match $plugin 'coordinateWeightUsesWindowCenter\s*=\s*true;' "The coordinate Gaussian center is not fixed to the window center."
Assert-Match $plugin '(?s)if\(!requested_normalize\)\s*\{\s*error_mesgs->report\("normalize_weights",.*?return false;' "Conflicting legacy normalization is not unconditionally rejected."
Assert-Match $plugin '(?s)if\(requested_zero_first\)\s*\{\s*error_mesgs->report\("zero_first_order_moment",.*?return false;' "Conflicting legacy first-moment centering is not rejected."
Assert-Match $plugin '(?s)if\(requested_invert\)\s*\{\s*error_mesgs->report\("invert_window_output",.*?return false;' "Conflicting legacy inversion is not rejected."
Assert-Match $plugin '(?s)if\(!center_mode\.empty\(\)\).*?normalized_center_mode != "window_center".*?error_mesgs->report\("weight_center_mode",.*?return false;' "Conflicting legacy center mode is not rejected."
if ($plugin -match 'parse_strict_bool_text\([^\r\n]*&(weightedOptions\.(normalize_weights|zero_first_order_moment)|invertWindowOutput)\)') {
    throw "A legacy boolean parameter still writes directly into fixed execution state."
}
foreach ($example in @('window_uniform.par', 'window_coordinate_gaussian.par', 'window_soft_high_porosity.par')) {
    [xml]$exampleXml = Get-Content -LiteralPath (Join-Path $RepoRoot "examples\$example") -Raw
    foreach ($fixedField in @('normalize_weights', 'zero_first_order_moment', 'invert_window_output', 'weight_center_mode')) {
        if ($null -ne $exampleXml.SelectSingleNode("/parameters/$fixedField")) {
            throw "The current example $example still serializes the fixed field $fixedField."
        }
    }
}
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
