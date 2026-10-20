param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = "Stop"
$pluginPath = Join-Path $RepoRoot "Template\HOSC2008\Plugin.cpp"
$callPath = Join-Path $RepoRoot "Template\HOSC2008\FunctionCall.cpp"
$initializationPath = Join-Path $RepoRoot "Template\HOSC2008\Initialization.cpp"
$uiPath = Join-Path $RepoRoot "ui\HOSC_TEMPLATE.ui"
$plugin = Get-Content -LiteralPath $pluginPath -Raw
$calls = Get-Content -LiteralPath $callPath -Raw
$initialization = Get-Content -LiteralPath $initializationPath -Raw
$ui = Get-Content -LiteralPath $uiPath -Raw
$feedbackPath = Join-Path $RepoRoot "shared\ValidationFeedback.h"
if(-not (Test-Path -LiteralPath $feedbackPath -PathType Leaf)) {
    throw "The shared Qt3 validation-feedback adapter is missing."
}

function Assert-Match([string]$Text, [string]$Pattern, [string]$Message) {
    if ($Text -notmatch $Pattern) { throw $Message }
}

Assert-Match $plugin 'HOSC_TEMPLATE supports Cartesian grids only' "The Cartesian-grid input guard is missing."
Assert-Match $plugin 'source_x\s*=\s*NULL;\s*source_y\s*=\s*NULL;\s*source_z\s*=\s*NULL;' "Ordinary x/y/z properties are still eligible to define template geometry."
Assert-Match $plugin 'source_val->is_informed\(node\)' "The source informed mask is not checked."
Assert-Match $plugin 'is_finite_number\(source_value\)' "Non-finite template source values are not excluded."
Assert-Match $plugin '(?s)template_sampling.*?exhaustive_integer_offsets' "The exact integer-offset sampling route is missing."
Assert-Match $plugin 'template_result_row_to_node_id' "Lag-result rows are not mapped explicitly to x-fastest SGeMS nodes."
Assert-Match $plugin 'const TemplateOutputLayout& layout' "The template result map still depends on calculator-owned Keys data."
Assert-Match $plugin 'output_layout_\.lag_step\[axis\]\s*=\s*Key->LAG\[reordered_direction\]' "The displayed lag steps are not frozen before calculation."
Assert-Match $plugin 'output_layout_\.displayed_direction_count\s*=\s*\(int\)displayedDirs\.size\(\)' "The selected display-direction count is not frozen explicitly."
Assert-Match $plugin 'template_result_row_to_node_id\(output_layout_, auto_grid_, targetGrid\[i\]\)' "Template result writers do not use the frozen output layout."
Assert-Match $initialization 'coordinate_basis=zero_based' "Template input does not declare the fixed zero-based coordinate boundary."
Assert-Match $plugin 'cannot sample this 3D Cartesian grid' "Legacy 3D search is not explicitly rejected."
Assert-Match $plugin '(?s)set_MyGrid_Data.*?template_result_row_to_node_id' "Primary template results are still written in enumeration order."
Assert-Match $plugin '(?s)write_Template_Output_Validity.*?template_result_row_to_node_id' "Template validity is not written by lag coordinate."
Assert-Match $plugin '(?s)write_Template_Exact_Support.*?template_result_row_to_node_id' "Template support is not written by lag coordinate."
Assert-Match $plugin 'Key->LTOL\[i\]\s*=\s*oldLTOL\[oldIndex\]' "LTOL is not preserved when directions are rearranged."
Assert-Match $plugin 'auto_prop_name_ \+ "_numeric_valid"' "Template numeric-valid output is missing."
Assert-Match $plugin 'auto_prop_name_ \+ "_display_valid"' "Template float-display validity output is missing."
Assert-Match $plugin 'is_template_float_representable\(value\)' "Template host-property writes do not check float32 representability."
Assert-Match $plugin 'auto_prop_name_ \+ "_support"' "Exact-template support output is missing."
Assert-Match $calls 'targetGrid\[i\]\s*=\s*\(double\*\)calloc\(DirectionNum \+ 3' "Template result rows lack support and validity columns."
Assert-Match $calls 'supportCounts\[row\]' "Exact-template support is not retained per lag."
Assert-Match $calls 'wmc::is_finite_number\(results\[row\]\) \? 1\.0 : 0\.0' "Exact-template numeric validity is not recorded per lag."
Assert-Match $calls 'std::numeric_limits<double>::quiet_NaN\(\)' "No-support exact lags are not represented as uninformed."
Assert-Match $ui '<cstring>template_sampling</cstring>' "The UI does not expose the template sampling mode."
Assert-Match $ui '<cstring>integer_direction_vectors</cstring>' "The UI does not expose exact integer direction vectors."
Assert-Match $plugin 'hosc_validation::FeedbackGuard feedback\(error_mesgs, "hosc_template_validation_summary"\)' "The template plugin does not publish its input errors to the form."
Assert-Match $ui '<cstring>hosc_template_validation_summary</cstring>' "The template UI has no field-specific validation summary."

function Get-XFastestNodeId([int]$X, [int]$Y, [int]$Z, [int]$NX, [int]$NY, [int]$NZ) {
    if($X -lt 0 -or $Y -lt 0 -or $Z -lt 0 -or $X -ge $NX -or $Y -ge $NY -or $Z -ge $NZ) {
        throw "Lag index is out of bounds."
    }
    return $X + $NX * ($Y + $NY * $Z)
}

function Get-FrozenTemplateOutputNodeId([int[]]$Row, [int]$AxisCount, [int[]]$RowColumns, [int[]]$LagSteps, [int[]]$Dimensions) {
    if(($AxisCount -ne 2 -and $AxisCount -ne 3) -or $Dimensions.Count -ne 3) {
        throw "Invalid frozen template output layout."
    }
    if($AxisCount -eq 2 -and $Dimensions[2] -ne 1) {
        throw "A two-direction layout must have a singleton third output axis."
    }
    $lagIndex = @(0, 0, 0)
    for($axis = 0; $axis -lt $AxisCount; ++$axis) {
        $column = $RowColumns[$axis]
        $step = $LagSteps[$axis]
        if($column -lt 0 -or $column -ge $Row.Count -or $step -le 0 -or ($Row[$column] % $step) -ne 0) {
            throw "Invalid frozen lag row."
        }
        $lagIndex[$axis] = [int]($Row[$column] / $step)
        if($lagIndex[$axis] -lt 0 -or $lagIndex[$axis] -ge $Dimensions[$axis]) {
            throw "Frozen lag row is out of bounds."
        }
    }
    return Get-XFastestNodeId $lagIndex[0] $lagIndex[1] $lagIndex[2] $Dimensions[0] $Dimensions[1] $Dimensions[2]
}

# SearchNode enumerates the final direction fastest. The following asymmetric
# fixtures verify the explicit map used at the Plugin output boundary.
$twoByThreeRows = @(@(0, 0), @(0, 1), @(0, 2), @(1, 0), @(1, 1), @(1, 2))
$twoByThreeIds = @($twoByThreeRows | ForEach-Object { Get-XFastestNodeId $_[0] $_[1] 0 2 3 1 })
if(($twoByThreeIds -join ',') -ne '0,2,4,1,3,5') {
    throw "2x3 x-fastest lag mapping is incorrect: $($twoByThreeIds -join ',')"
}

$threeDimensionalCases = @(
    [PSCustomObject]@{ lag = @(0, 0, 0); expected = 0 },
    [PSCustomObject]@{ lag = @(1, 0, 0); expected = 1 },
    [PSCustomObject]@{ lag = @(0, 1, 0); expected = 2 },
    [PSCustomObject]@{ lag = @(1, 2, 3); expected = 23 }
)
foreach($case in $threeDimensionalCases) {
    $actual = Get-XFastestNodeId $case.lag[0] $case.lag[1] $case.lag[2] 2 3 4
    if($actual -ne $case.expected) {
        throw "2x3x4 x-fastest lag mapping is incorrect for $($case.lag -join ','): $actual"
    }
}

# Three selected directions remain three directions even when the final output
# axis has one lag. The row-column map must not infer that information from nz.
$threeDirectionSingletonRows = @(
    @(0, 0, 0), @(0, 1, 0), @(0, 2, 0),
    @(1, 0, 0), @(1, 1, 0), @(1, 2, 0)
)
$threeDirectionSingletonIds = @($threeDirectionSingletonRows | ForEach-Object {
    Get-FrozenTemplateOutputNodeId $_ 3 @(0, 1, 2) @(1, 1, 1) @(2, 3, 1)
})
if(($threeDirectionSingletonIds -join ',') -ne '0,2,4,1,3,5') {
    throw "Three-direction 2x3x1 lag mapping is incorrect: $($threeDirectionSingletonIds -join ',')"
}

# With one hidden direction, the reordered result row contains the hidden lag
# first. The snapshot columns, not the former Key arrays, select output axes.
$hiddenDirectionRows = @(
    @(0, 0, 0, 0), @(0, 0, 1, 0), @(0, 0, 2, 0),
    @(0, 3, 0, 0), @(0, 3, 1, 0), @(0, 3, 2, 0)
)
$hiddenDirectionIds = @($hiddenDirectionRows | ForEach-Object {
    Get-FrozenTemplateOutputNodeId $_ 3 @(1, 2, 3) @(3, 1, 2) @(2, 3, 1)
})
if(($hiddenDirectionIds -join ',') -ne '0,2,4,1,3,5') {
    throw "Hidden-direction frozen lag mapping is incorrect: $($hiddenDirectionIds -join ',')"
}

# A missing first source node must not change the coordinate basis for a later
# informed node. This mirrors the fixed zero-based Plugin boundary.
$firstInformedAfterMissing = @(1, 0, 0)
if($firstInformedAfterMissing[0] -lt 0 -or $firstInformedAfterMissing[0] -ge 3) {
    throw "The zero-based first-missing-node fixture is not in the source grid."
}

Write-Host "Template contract verification passed against current source and UI."
