#Requires -Version 5.1
param(
    [Parameter(Mandatory = $true)][string]$SGeMSRoot,
    [string]$FixturePath = '',
    [string]$ExpectedWindowPath = '',
    [string]$ExampleDirectory = '',
    [string]$EvidenceRoot = '',
    [int]$TimeoutSeconds = 120,
    [switch]$PrepareOnly
)

$ErrorActionPreference = 'Stop'
if (!$FixturePath) { $FixturePath = Join-Path $PSScriptRoot '..\tests\fixtures\manual\selftest_3d_cartesian.sgems' }
if (!$ExpectedWindowPath) { $ExpectedWindowPath = Join-Path $PSScriptRoot '..\tests\fixtures\manual\expected_window.csv' }
if (!$ExampleDirectory) { $ExampleDirectory = Join-Path $PSScriptRoot '..\examples' }
$invariant = [Globalization.CultureInfo]::InvariantCulture
if ($TimeoutSeconds -le 0) { throw 'TimeoutSeconds must be positive.' }
$sgemsExe = [IO.Path]::GetFullPath((Join-Path $SGeMSRoot 'sgems.exe'))
foreach ($required in @($sgemsExe, $FixturePath, $ExpectedWindowPath)) {
    if (!(Test-Path -LiteralPath $required -PathType Leaf)) { throw "Missing required file: $required" }
}
if (!(Test-Path -LiteralPath $ExampleDirectory -PathType Container)) { throw "Missing GUI parameter-file directory: $ExampleDirectory" }
# Other isolated copies may coexist. Never close or reuse an existing process.
foreach ($running in @(Get-Process sgems -ErrorAction SilentlyContinue)) {
    try { $runningPath = $running.Path } catch { throw "Cannot identify SGeMS PID $($running.Id); refusing to start." }
    if ($runningPath -and [string]::Equals([IO.Path]::GetFullPath($runningPath), $sgemsExe, [StringComparison]::OrdinalIgnoreCase)) {
        throw "SGeMS is already running from this root (PID $($running.Id)). Supply a separate copied host root."
    }
}
if (!$EvidenceRoot) {
    $EvidenceRoot = Join-Path ([IO.Path]::GetTempPath()) ('hosc_validation_' + [Guid]::NewGuid().ToString('N'))
}
$EvidenceRoot = [IO.Path]::GetFullPath($EvidenceRoot)
if (Test-Path -LiteralPath $EvidenceRoot) { throw "Evidence destination already exists: $EvidenceRoot" }
New-Item -ItemType Directory -Path $EvidenceRoot | Out-Null
$auditRoot = Join-Path $EvidenceRoot 'audit'
New-Item -ItemType Directory -Path $auditRoot | Out-Null

function Read-BigEndianUInt32([IO.BinaryReader]$Reader) {
    $bytes = $Reader.ReadBytes(4)
    if ($bytes.Length -ne 4) { throw 'Unexpected end of S-GeMS file.' }
    [Array]::Reverse($bytes)
    return [BitConverter]::ToUInt32($bytes, 0)
}
function Read-BigEndianSingle([IO.BinaryReader]$Reader) {
    $bytes = $Reader.ReadBytes(4)
    if ($bytes.Length -ne 4) { throw 'Unexpected end of S-GeMS property data.' }
    [Array]::Reverse($bytes)
    return [double][BitConverter]::ToSingle($bytes, 0)
}
function Read-CString([IO.BinaryReader]$Reader) {
    $length = Read-BigEndianUInt32 $Reader
    if ($length -lt 1 -or $length -gt 1048576) { throw "Invalid S-GeMS string length: $length" }
    $bytes = $Reader.ReadBytes($length)
    if ($bytes.Length -ne $length -or $bytes[$length - 1] -ne 0) { throw 'Invalid S-GeMS string terminator.' }
    return [Text.Encoding]::ASCII.GetString($bytes, 0, $length - 1)
}
function Read-CartesianGrid([string]$Path) {
    $stream = [IO.File]::OpenRead($Path)
    try {
        $reader = New-Object IO.BinaryReader($stream)
        if ((Read-BigEndianUInt32 $reader) -ne [uint32]::Parse('B211175D', [Globalization.NumberStyles]::HexNumber)) { throw "Invalid grid magic: $Path" }
        if ((Read-CString $reader) -ne 'Cgrid') { throw "Expected Cartesian grid: $Path" }
        $name = Read-CString $reader
        if ((Read-BigEndianUInt32 $reader) -ne 100) { throw 'Unsupported grid file version.' }
        $nx = Read-BigEndianUInt32 $reader
        $ny = Read-BigEndianUInt32 $reader
        $nz = Read-BigEndianUInt32 $reader
        if ($nx -lt 1 -or $ny -lt 1 -or $nz -lt 1 -or ($nx * $ny * $nz) -gt 10000) { throw 'Unexpected test-grid dimensions.' }
        $geometry = @()
        for ($i = 0; $i -lt 6; ++$i) { $geometry += Read-BigEndianSingle $reader }
        $count = Read-BigEndianUInt32 $reader
        if ($count -gt 100) { throw 'Unexpected test-grid property count.' }
        $names = @()
        for ($i = 0; $i -lt $count; ++$i) { $names += Read-CString $reader }
        $properties = @{}
        foreach ($property in $names) {
            $values = New-Object 'System.Collections.Generic.List[double]'
            for ($i = 0; $i -lt ($nx * $ny * $nz); ++$i) { $values.Add((Read-BigEndianSingle $reader)) }
            $properties[$property] = $values.ToArray()
        }
        if ($stream.Position -ne $stream.Length) { throw "Unexpected trailing grid data: $Path" }
        return [PSCustomObject]@{ Name = $name; NX = $nx; NY = $ny; NZ = $nz; Geometry = $geometry; Properties = $properties }
    }
    finally { $stream.Dispose() }
}
function New-ParameterXml([string]$Algorithm, [Collections.IDictionary]$Values) {
    $parts = New-Object 'System.Collections.Generic.List[string]'
    $parts.Add('<parameters><algorithm name="' + $Algorithm + '" />')
    foreach ($key in $Values.Keys) {
        $parts.Add('<' + $key + ' value="' + [Security.SecurityElement]::Escape([string]$Values[$key]) + '" />')
    }
    $parts.Add('</parameters>')
    return $parts -join ''
}
function Assert-Close([double]$Actual, [double]$Expected, [double]$Tolerance, [string]$Label) {
    if ([double]::IsNaN($Actual) -or [double]::IsInfinity($Actual) -or [Math]::Abs($Actual - $Expected) -gt $Tolerance) {
        throw "$Label expected $Expected, got $Actual (absolute tolerance $Tolerance)."
    }
}
function Read-KeyValues([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $separator = $line.IndexOf('=')
        if ($separator -gt 0) { $values[$line.Substring(0, $separator)] = $line.Substring($separator + 1) }
    }
    return $values
}

$source = Read-CartesianGrid $FixturePath
if ($source.Name -ne 'test' -or $source.NX -ne 5 -or $source.NY -ne 5 -or $source.NZ -ne 3 -or !$source.Properties.ContainsKey('phi')) {
    throw 'This regression requires the test / phi Cartesian fixture with dimensions 5 x 5 x 3.'
}
$template = [ordered]@{
    Target_Button = 'Cumulant'; Algorithm_Mode = 'AVL'; Source_Grid = 'test'; Source_Property = 'phi'; props = 'invalid_template';
    xmn_ = '0'; ymn_ = '0'; zmn_ = '0';
    NOD = '1,1,1'; NLAG = '6,12,12'; LAG = '4,4,4'; LTOL = '0,0,0';
    AZM = '180,90,0'; ATOL = '0,0,0'; BANDWH = '0,0,0'; DIP = '0,0,0'; DTOL = '0,0,0'; BANDWD = '0,0,0';
    fst_dir = '3'; snd_dir = '2'; trd_dir = '0'; cross_pos = '1'; Output_3D_Slice_Enable = 'Disable';
    Weighted_Enable = 'Disable'; Weight_Mode = 'uniform'; Normalize_Weights = 'true'; Zero_First_Order_Moment = 'false';
    Weight_Bandwidth = '1'; Weight_Minimum = '0'; Weight_Reference_Value = '0'; Weight_Center_X = '0'; Weight_Center_Y = '0'; Weight_Center_Z = '0';
    template_sampling = 'legacy_search'; integer_direction_vectors = '1,0;0,1'
}
$invalidTemplateXml = New-ParameterXml 'HOSC_TEMPLATE' $template
$template['props'] = 'valid_template'
$template['NLAG'] = '3,3,3'; $template['LAG'] = '1,1,1'
$template['fst_dir'] = '1'; $template['snd_dir'] = '2'; $template['trd_dir'] = '3'
$template['Weighted_Enable'] = 'Disable'
$template['template_sampling'] = 'exhaustive_integer_offsets'
$template['integer_direction_vectors'] = '1,0,0;0,1,0;0,0,1'
$validTemplateXml = New-ParameterXml 'HOSC_TEMPLATE' $template
$window = [ordered]@{
    Target_Button = 'Cumulant'; Algorithm_Button = 'Bell'; Source_Grid = 'test'; Source_Property = 'phi'; props = 'invalid_window';
    moment_order = '4'; windowl = '0'; windoww = '3'; windowh = '3'; stepx = '1'; stepy = '1'; stepz = '1';
    xmn_ = '0'; ymn_ = '0'; zmn_ = '0'; invert_window_output = 'false';
    enable_weighted_moment = 'false'; weight_mode = 'uniform'; normalize_weights = 'true';
    center_cumulant_input_for_order_ge_2 = 'true'; enable_percentile_clipping = 'false'; percentile_vmin = '2'; percentile_vmax = '98';
    output_standardized_cumulant = 'true'; variance_epsilon = '1e-12'; output_diagnostics = 'true';
    write_audit_sidecar = 'true'; audit_output_directory = $auditRoot
}
$invalidWindowXml = New-ParameterXml 'HOSC_WINDOW' $window
$window['windowl'] = '3'; $window['props'] = 'valid_window'
$validWindowXml = New-ParameterXml 'HOSC_WINDOW' $window
# Exercise the removed fields as legacy PAR input. Each invalid request must
# fail before creating output, then the normal calls below must still succeed.
$policyCases = @()
foreach ($spec in @(
    @{Name='normalize_enabled'; Field='normalize_weights'; Value='false'; Enabled='true'; Mode='uniform'},
    @{Name='normalize_disabled'; Field='normalize_weights'; Value='false'; Enabled='false'; Mode='uniform'},
    @{Name='legacy_centering'; Field='zero_first_order_moment'; Value='true'; Enabled='true'; Mode='uniform'},
    @{Name='inverted_statistic'; Field='invert_window_output'; Value='true'; Enabled='false'; Mode='uniform'},
    @{Name='parameter_center'; Field='weight_center_mode'; Value='parameter'; Enabled='true'; Mode='coordinate_gaussian'},
    @{Name='parameter_center_uniform'; Field='weight_center_mode'; Value='parameter'; Enabled='false'; Mode='uniform'}
)) {
    $parameters = [ordered]@{}
    foreach ($key in $window.Keys) { $parameters[$key] = $window[$key] }
    $parameters['props'] = 'rejected_' + $spec.Name
    $parameters['enable_weighted_moment'] = $spec.Enabled
    $parameters['weight_mode'] = $spec.Mode
    $parameters[$spec.Field] = $spec.Value
    $policyCases += [PSCustomObject]@{Name=$spec.Name; Algorithm='HOSC_WINDOW'; Field=$spec.Field; Xml=(New-ParameterXml 'HOSC_WINDOW' $parameters)}
}
# A real one-layer Cartesian fixture covers the former legacy-2D bypass.
function Write-BigEndianNumber($Writer, $Value, [switch]$Single) {
    if ($Single) { $bytes = [BitConverter]::GetBytes([single]$Value) }
    else { $bytes = [BitConverter]::GetBytes([uint32]$Value) }
    [Array]::Reverse($bytes); $Writer.Write($bytes)
}
function Write-GridString($Writer, [string]$Value) {
    $bytes = [Text.Encoding]::ASCII.GetBytes($Value + [char]0)
    Write-BigEndianNumber $Writer $bytes.Length
    $Writer.Write($bytes)
}
$fixture2d = Join-Path $EvidenceRoot 'input_2d.sgems'
$fixtureWriter = New-Object IO.BinaryWriter([IO.File]::Create($fixture2d))
try {
    Write-BigEndianNumber $fixtureWriter ([uint32]::Parse('B211175D',[Globalization.NumberStyles]::HexNumber))
    Write-GridString $fixtureWriter 'Cgrid'; Write-GridString $fixtureWriter 'test2d'
    foreach ($number in @(100,5,5,1)) { Write-BigEndianNumber $fixtureWriter $number }
    foreach ($number in $source.Geometry) { Write-BigEndianNumber $fixtureWriter $number -Single }
    Write-BigEndianNumber $fixtureWriter 1; Write-GridString $fixtureWriter 'phi'
    foreach ($number in $source.Properties['phi'][0..24]) { Write-BigEndianNumber $fixtureWriter $number -Single }
} finally { $fixtureWriter.Dispose() }
if ((Read-CartesianGrid $fixture2d).NZ -ne 1) { throw 'Invalid 2D fixture.' }
foreach ($mode in @('legacy_search','exhaustive_integer_offsets')) {
    $parameters = [ordered]@{}
    foreach ($key in $template.Keys) { $parameters[$key] = $template[$key] }
    $parameters['Weighted_Enable'] = 'Enable'; $parameters['props'] = 'rejected_' + $mode
    $parameters['template_sampling'] = $mode
    if ($mode -eq 'legacy_search') { $parameters['Source_Grid'] = 'test2d' }
    $policyCases += [PSCustomObject]@{Name=$mode; Algorithm='HOSC_TEMPLATE'; Field='Weighted_Enable'; Xml=(New-ParameterXml 'HOSC_TEMPLATE' $parameters)}
}
$templateExport = Join-Path $EvidenceRoot 'valid_template.sgems'
$windowExport = Join-Path $EvidenceRoot 'valid_window.sgems'
$commandPath = Join-Path $EvidenceRoot 'validation_sequence.cmd'
$stdoutPath = Join-Path $EvidenceRoot 'host.stdout.log'
$stderrPath = Join-Path $EvidenceRoot 'host.stderr.log'
$commands = @(
    "LoadObjectFromFile $FixturePath`::s-gems",
    "LoadObjectFromFile $fixture2d`::s-gems"
)
foreach ($case in $policyCases) {
    $commands += 'RunGeostatAlgorithm ' + $case.Algorithm + '::/GeostatParamUtils/XML::' + $case.Xml
}
$commands += @(
    "RunGeostatAlgorithm HOSC_TEMPLATE::/GeostatParamUtils/XML::$invalidTemplateXml",
    "RunGeostatAlgorithm HOSC_TEMPLATE::/GeostatParamUtils/XML::$validTemplateXml",
    "SaveGeostatGrid myGrid::$templateExport`::s-gems",
    "RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::$invalidWindowXml",
    "RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::$validWindowXml",
    "SaveGeostatGrid myWM::$windowExport`::s-gems"
)
# Execute the actual files shipped for GUI loading. Do not synthesize their
# parameter names, enable an audit sidecar, or substitute numerical settings.
# In particular this catches a GUI PAR using the wrong algorithm selector key.
$executedParameters = Join-Path $EvidenceRoot 'executed_parameters'
New-Item -ItemType Directory -Path $executedParameters | Out-Null
$exampleSpecs = @(
    # The successful recovery calls above own myGrid and myWM. Each following
    # call creates the next unique grid; export that grid, not the first result.
    [PSCustomObject]@{ Name = 'template_3d.par'; Algorithm = 'HOSC_TEMPLATE'; CalculatorKey = 'Algorithm_Mode'; Calculator = 'AVL'; WeightMode = 'uniform'; Order = 4; Grid = 'myGrid_1' },
    [PSCustomObject]@{ Name = 'window_uniform.par'; Algorithm = 'HOSC_WINDOW'; CalculatorKey = 'Algorithm_Button'; Calculator = 'Bell'; WeightMode = 'uniform'; Order = 4; Grid = 'myWM_1' },
    [PSCustomObject]@{ Name = 'window_coordinate_gaussian.par'; Algorithm = 'HOSC_WINDOW'; CalculatorKey = 'Algorithm_Button'; Calculator = 'Bell'; WeightMode = 'coordinate_gaussian'; Order = 6; Grid = 'myWM_2' },
    [PSCustomObject]@{ Name = 'window_soft_high_porosity.par'; Algorithm = 'HOSC_WINDOW'; CalculatorKey = 'Algorithm_Button'; Calculator = 'Bell'; WeightMode = 'soft_high_porosity'; Order = 6; Grid = 'myWM_3' }
)
$exampleCases = @()
$usedOutputNames = @('valid_template', 'valid_window', 'invalid_template', 'invalid_window')
foreach ($spec in $exampleSpecs) {
    $parameterPath = [IO.Path]::GetFullPath((Join-Path $ExampleDirectory $spec.Name))
    if (!(Test-Path -LiteralPath $parameterPath -PathType Leaf)) { throw "Missing GUI parameter example: $parameterPath" }
    $document = New-Object Xml.XmlDocument
    $document.XmlResolver = $null
    $document.PreserveWhitespace = $false
    $document.Load($parameterPath)
    if ($document.DocumentType) { throw "Document types are not allowed in parameter examples: $($spec.Name)" }
    $algorithmNode = $document.SelectSingleNode('/parameters/algorithm')
    $calculatorNode = $document.SelectSingleNode('/parameters/' + $spec.CalculatorKey)
    $propertyNode = $document.SelectSingleNode('/parameters/props')
    if (!$algorithmNode -or $algorithmNode.GetAttribute('name') -ne $spec.Algorithm) { throw "Unexpected algorithm in $($spec.Name)" }
    if (!$calculatorNode -or $calculatorNode.GetAttribute('value') -ne $spec.Calculator) {
        throw "$($spec.Name) must provide $($spec.CalculatorKey)=$($spec.Calculator). The actual GUI selector key is required."
    }
    if (!$propertyNode -or [string]::IsNullOrWhiteSpace($propertyNode.GetAttribute('value'))) { throw "Missing output property in $($spec.Name)" }
    $propertyName = $propertyNode.GetAttribute('value')
    if ($usedOutputNames -contains $propertyName) { throw "Duplicate example output property: $propertyName" }
    $usedOutputNames += $propertyName
    if ($spec.Algorithm -eq 'HOSC_WINDOW') {
        $orderNode = $document.SelectSingleNode('/parameters/moment_order')
        $weightNode = $document.SelectSingleNode('/parameters/weight_mode')
        if (!$orderNode -or $orderNode.GetAttribute('value') -ne [string]$spec.Order -or !$weightNode -or $weightNode.GetAttribute('value') -ne $spec.WeightMode) {
            throw "Unexpected weight mode or order in $($spec.Name); expected $($spec.WeightMode), order $($spec.Order)."
        }
    }
    $baseName = [IO.Path]::GetFileNameWithoutExtension($spec.Name)
    $exportPath = Join-Path $EvidenceRoot ($baseName + '.sgems')
    Copy-Item -LiteralPath $parameterPath -Destination (Join-Path $executedParameters $spec.Name)
    $xml = $document.OuterXml
    if ($xml -match '[\r\n]') { throw "Parameter XML cannot be represented on one command line: $($spec.Name)" }
    $commands += 'RunGeostatAlgorithm ' + $spec.Algorithm + '::/GeostatParamUtils/XML::' + $xml
    $commands += 'SaveGeostatGrid ' + $spec.Grid + '::' + $exportPath + '::s-gems'
    $exampleCases += [PSCustomObject]@{
        Name = $spec.Name; Algorithm = $spec.Algorithm; CalculatorKey = $spec.CalculatorKey; Calculator = $spec.Calculator;
        WeightMode = $spec.WeightMode; Order = $spec.Order; Grid = $spec.Grid; Property = $propertyName; ExportPath = $exportPath;
        ParameterSha256 = (Get-FileHash -LiteralPath $parameterPath -Algorithm SHA256).Hash.ToLowerInvariant()
    }
}
if (($commands -join "`n") -match '[^\x00-\x7F]') { throw 'SGeMS command-file mode requires ASCII fixture and evidence paths for this test.' }
[IO.File]::WriteAllLines($commandPath, $commands, [Text.Encoding]::ASCII)
if ($PrepareOnly) {
    Write-Host "Prepared validation regression without starting SGeMS: $EvidenceRoot"
    return
}

$startInfo = New-Object Diagnostics.ProcessStartInfo
$startInfo.FileName = $sgemsExe
$startInfo.Arguments = '"{0}"' -f $commandPath
$startInfo.WorkingDirectory = $SGeMSRoot
$startInfo.UseShellExecute = $false
$startInfo.CreateNoWindow = $true
$startInfo.WindowStyle = [Diagnostics.ProcessWindowStyle]::Hidden
$startInfo.RedirectStandardOutput = $true
$startInfo.RedirectStandardError = $true
$process = New-Object Diagnostics.Process
$process.StartInfo = $startInfo
$stdoutTask = $null
$stderrTask = $null
try {
    # Keep the original process handle. Start-Process -PassThru in Windows
    # PowerShell 5.1 can lose ExitCode for a fast-exiting child GUI executable.
    if (!$process.Start()) { throw 'Could not start the isolated SGeMS process.' }
    $processHandle = $process.Handle
    $hostPid = $process.Id
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $deadline = [DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
    while (!$process.WaitForExit(1000)) {
        if ([DateTime]::UtcNow -ge $deadline) {
            # This object is the child created above, never a user-owned GUI.
            $process.Kill()
            $process.WaitForExit()
            throw "Regression timed out (possible modal dialog). Only child PID $($process.Id) was stopped. Evidence: $EvidenceRoot"
        }
    }
    $process.WaitForExit()
    $exitCode = $process.ExitCode
    if ($null -eq $exitCode) { throw "The isolated SGeMS exit code could not be obtained. Evidence: $EvidenceRoot" }
    if ($exitCode -ne 0) { throw "Child host exited with $exitCode. Evidence: $EvidenceRoot" }
}
finally {
    if ($stdoutTask) { [IO.File]::WriteAllText($stdoutPath, $stdoutTask.GetAwaiter().GetResult(), (New-Object Text.UTF8Encoding($false))) }
    if ($stderrTask) { [IO.File]::WriteAllText($stderrPath, $stderrTask.GetAwaiter().GetResult(), (New-Object Text.UTF8Encoding($false))) }
    $process.Dispose()
}

$log = (Get-Content -LiteralPath $stdoutPath -Raw) + "`n" + (Get-Content -LiteralPath $stderrPath -Raw)
$issues = New-Object 'System.Collections.Generic.List[object]'
foreach ($match in [regex]::Matches($log, '<(?<field>[A-Za-z_][A-Za-z0-9_]*)\s+diagnostic="(?<message>[^\r\n]*)"\s*/>')) {
    $issues.Add([PSCustomObject]@{ Field = $match.Groups['field'].Value; Message = [Net.WebUtility]::HtmlDecode($match.Groups['message'].Value) })
}
foreach ($field in @('template_sampling', 'integer_direction_vectors', 'windowl', 'Weighted_Enable', 'normalize_weights', 'zero_first_order_moment', 'invert_window_output', 'weight_center_mode')) {
    $expectedCount = @($policyCases | Where-Object { $_.Field -eq $field }).Count
    if ($field -in @('template_sampling', 'integer_direction_vectors', 'windowl')) { $expectedCount++ }
    if (@($issues | Where-Object { $_.Field -eq $field }).Count -ne $expectedCount) { throw "Expected $expectedCount diagnostic(s) for $field. See $stderrPath" }
}
if ($issues.Count -ne (3 + $policyCases.Count)) { throw "Unexpected total field diagnostic count: $($issues.Count)." }
if ([regex]::Matches($log, 'Command RunGeostatAlgorithm could not be performed').Count -ne (2 + $policyCases.Count)) { throw 'Only the deliberately invalid algorithm calls should fail.' }
$issues | Export-Csv -LiteralPath (Join-Path $EvidenceRoot 'reported_issues.csv') -NoTypeInformation -Encoding UTF8

$output = Read-CartesianGrid $templateExport
if ($output.Name -ne 'myGrid' -or $output.NX -ne 3 -or $output.NY -ne 3 -or $output.NZ -ne 3) { throw 'The valid exact-template result is not the expected 3 x 3 x 3 lag grid.' }
foreach ($property in @('valid_template', 'valid_template_support', 'valid_template_numeric_valid', 'valid_template_display_valid')) {
    if (!$output.Properties.ContainsKey($property)) { throw "Missing exact-template property: $property" }
}
if ($output.Properties.ContainsKey('invalid_template')) { throw 'The invalid template unexpectedly created an output property.' }
$templateChecks = New-Object 'System.Collections.Generic.List[object]'
for ($lz = 0; $lz -lt 3; ++$lz) {
    for ($ly = 0; $ly -lt 3; ++$ly) {
        for ($lx = 0; $lx -lt 3; ++$lx) {
            $node = $lx + 3 * ($ly + 3 * $lz)
            $samples = New-Object 'System.Collections.Generic.List[object]'
            [double[]]$means = @(0, 0, 0, 0)
            for ($z = 0; $z -lt 3 - $lz; ++$z) {
                for ($y = 0; $y -lt 5 - $ly; ++$y) {
                    for ($x = 0; $x -lt 5 - $lx; ++$x) {
                        $anchor = $x + 5 * ($y + 5 * $z)
                        [double[]]$sample = @($source.Properties['phi'][$anchor], $source.Properties['phi'][$anchor + $lx], $source.Properties['phi'][$anchor + 5 * $ly], $source.Properties['phi'][$anchor + 25 * $lz])
                        $samples.Add($sample)
                        for ($i = 0; $i -lt 4; ++$i) { $means[$i] += $sample[$i] }
                    }
                }
            }
            for ($i = 0; $i -lt 4; ++$i) { $means[$i] /= $samples.Count }
            $m4 = 0.0; [double[]]$cov = @(0, 0, 0, 0, 0, 0)
            foreach ($sample in $samples) {
                [double[]]$v = @(0, 0, 0, 0)
                for ($i = 0; $i -lt 4; ++$i) { $v[$i] = $sample[$i] - $means[$i] }
                $m4 += $v[0] * $v[1] * $v[2] * $v[3]
                $cov[0] += $v[0] * $v[1]; $cov[1] += $v[2] * $v[3]
                $cov[2] += $v[0] * $v[2]; $cov[3] += $v[1] * $v[3]
                $cov[4] += $v[0] * $v[3]; $cov[5] += $v[1] * $v[2]
            }
            $expected = $m4 / $samples.Count - ($cov[0] * $cov[1] + $cov[2] * $cov[3] + $cov[4] * $cov[5]) / ($samples.Count * $samples.Count)
            Assert-Close $output.Properties['valid_template_support'][$node] $samples.Count 0 "template support node $node"
            Assert-Close $output.Properties['valid_template_numeric_valid'][$node] 1 0 "template numeric validity node $node"
            Assert-Close $output.Properties['valid_template_display_valid'][$node] 1 0 "template display validity node $node"
            Assert-Close $output.Properties['valid_template'][$node] $expected 1e-9 "template float output node $node"
            $templateChecks.Add([PSCustomObject]@{ node_id = $node; lag_x = $lx; lag_y = $ly; lag_z = $lz; support = $samples.Count; expected = $expected; actual = $output.Properties['valid_template'][$node] })
        }
    }
}
$templateChecks | Export-Csv -LiteralPath (Join-Path $EvidenceRoot 'template_result_checks.csv') -NoTypeInformation -Encoding UTF8

$sidecars = @(Get-ChildItem -LiteralPath $auditRoot -Directory | Where-Object { $_.Name -like 'HOSC_WINDOW_*' })
if ($sidecars.Count -ne 1) { throw "Expected one successful window audit directory; got $($sidecars.Count)." }
$status = Read-KeyValues (Join-Path $sidecars[0].FullName 'run_status.txt')
$config = Read-KeyValues (Join-Path $sidecars[0].FullName 'resolved_config.txt')
if ($status['state'] -ne 'completed' -or $config['output_property'] -ne 'valid_window' -or [int]$config['process_id'] -ne $hostPid) { throw 'Window did not recover successfully in the same host process.' }
$rows = @(Import-Csv -LiteralPath (Join-Path $sidecars[0].FullName 'window_results.csv'))
$reference = @{}
foreach ($row in Import-Csv -LiteralPath $ExpectedWindowPath | Where-Object { $_.property -eq 'phi' -and $_.weight_mode -eq 'uniform' -and $_.order -eq '4' }) { $reference[[int]$row.node_id] = $row }
if ($rows.Count -ne 9 -or $reference.Count -ne 9) { throw 'Window or reference does not contain nine expected rows.' }
foreach ($row in $rows) {
    $node = [int]$row.node_id
    if (!$reference.ContainsKey($node)) { throw "Unexpected window node: $node" }
    if ($row.numeric_valid -ne '1' -or $row.display_valid -ne '1' -or $row.lambda_valid -ne '1') { throw "Invalid window result at node $node" }
    Assert-Close ([double]::Parse($row.raw_value, $invariant)) ([double]::Parse($reference[$node].raw_value, $invariant)) 1e-12 "window cumulant node $node"
    Assert-Close ([double]::Parse($row.lambda, $invariant)) ([double]::Parse($reference[$node].lambda, $invariant)) 1e-6 "window lambda node $node"
    Assert-Close ([double]$row.support) 27 0 "window support node $node"
}
$windowGrid = Read-CartesianGrid $windowExport
if ($windowGrid.NX -ne 3 -or $windowGrid.NY -ne 3 -or $windowGrid.NZ -ne 1 -or !$windowGrid.Properties.ContainsKey('valid_window')) {
    throw 'The exported window grid does not contain the expected 3 x 3 x 1 result.'
}
foreach ($row in $rows) {
    $node = [int]$row.node_id
    Assert-Close $windowGrid.Properties['valid_window'][$node] ([double]::Parse($row.raw_value, $invariant)) 1e-9 "exported window value node $node"
}
$exampleChecks = New-Object 'System.Collections.Generic.List[object]'
$exampleResults = @()
$allExpectedWindowRows = @(Import-Csv -LiteralPath $ExpectedWindowPath)
foreach ($case in $exampleCases) {
    $grid = Read-CartesianGrid $case.ExportPath
    $prefix = $case.Property
    $maxRawError = 0.0; $maxLambdaError = 0.0; $maxVarianceError = 0.0; $maxEssError = 0.0
    if ($case.Algorithm -eq 'HOSC_TEMPLATE') {
        if ($grid.Name -ne $case.Grid -or $grid.NX -ne 3 -or $grid.NY -ne 3 -or $grid.NZ -ne 3) { throw "Unexpected lag grid for $($case.Name)" }
        foreach ($suffix in @('', '_support', '_numeric_valid', '_display_valid')) {
            if (!$grid.Properties.ContainsKey($prefix + $suffix)) { throw "Missing property $prefix$suffix in $($case.Name)" }
        }
        foreach ($expectedRow in $templateChecks) {
            $node = [int]$expectedRow.node_id
            $actual = $grid.Properties[$prefix][$node]
            Assert-Close $actual $expectedRow.expected 1e-9 "$($case.Name) cumulant node $node"
            Assert-Close $grid.Properties[$prefix + '_support'][$node] $expectedRow.support 0 "$($case.Name) support node $node"
            Assert-Close $grid.Properties[$prefix + '_numeric_valid'][$node] 1 0 "$($case.Name) numeric validity node $node"
            Assert-Close $grid.Properties[$prefix + '_display_valid'][$node] 1 0 "$($case.Name) display validity node $node"
            $errorValue = [Math]::Abs($actual - $expectedRow.expected)
            $maxRawError = [Math]::Max($maxRawError, $errorValue)
            $exampleChecks.Add([PSCustomObject]@{ parameter_file = $case.Name; property = $prefix; weight_mode = $case.WeightMode; order = $case.Order; node_id = $node; expected_raw = $expectedRow.expected; actual_raw = $actual; expected_kappa2 = ''; actual_kappa2 = ''; expected_lambda = ''; actual_lambda = ''; expected_ess = ''; actual_ess = ''; support = $expectedRow.support })
        }
        $checkedNodes = 27
        $referenceDescription = 'Independent centered fourth joint moment minus the three covariance products'
    }
    else {
        if ($grid.Name -ne $case.Grid -or $grid.NX -ne 3 -or $grid.NY -ne 3 -or $grid.NZ -ne 1) { throw "Unexpected window grid for $($case.Name)" }
        foreach ($suffix in @('', '_kappa2', '_lambda', '_support', '_ess', '_numeric_valid', '_display_valid', '_lambda_valid')) {
            if (!$grid.Properties.ContainsKey($prefix + $suffix)) { throw "Missing property $prefix$suffix in $($case.Name)" }
        }
        $expectedRows = @($allExpectedWindowRows | Where-Object { $_.property -eq 'phi' -and $_.weight_mode -eq $case.WeightMode -and [int]$_.order -eq $case.Order })
        if ($expectedRows.Count -ne 9 -or @($expectedRows | Select-Object -ExpandProperty node_id -Unique).Count -ne 9) { throw "Expected nine distinct reference nodes for $($case.Name)" }
        foreach ($expectedRow in $expectedRows) {
            $node = [int]$expectedRow.node_id
            if ($node -lt 0 -or $node -ge 9) { throw "Unexpected reference node $node for $($case.Name)" }
            $expectedRaw = [double]::Parse($expectedRow.raw_value, $invariant)
            $expectedVariance = [double]::Parse($expectedRow.kappa2, $invariant)
            $expectedLambda = [double]::Parse($expectedRow.lambda, $invariant)
            $expectedEss = [double]::Parse($expectedRow.ess, $invariant)
            $actualRaw = $grid.Properties[$prefix][$node]
            $actualVariance = $grid.Properties[$prefix + '_kappa2'][$node]
            $actualLambda = $grid.Properties[$prefix + '_lambda'][$node]
            $actualEss = $grid.Properties[$prefix + '_ess'][$node]
            # SGeMS exports single-precision properties; the reference CSV uses
            # double precision. Allow float rounding, not changes of statistic.
            Assert-Close $actualRaw $expectedRaw ([Math]::Max(1e-12, 2e-7 * [Math]::Abs($expectedRaw))) "$($case.Name) cumulant node $node"
            Assert-Close $actualVariance $expectedVariance ([Math]::Max(1e-10, 2e-7 * [Math]::Abs($expectedVariance))) "$($case.Name) variance node $node"
            Assert-Close $actualLambda $expectedLambda ([Math]::Max(1e-6, 2e-7 * [Math]::Abs($expectedLambda))) "$($case.Name) lambda node $node"
            Assert-Close $actualEss $expectedEss 1e-5 "$($case.Name) ESS node $node"
            Assert-Close $grid.Properties[$prefix + '_support'][$node] ([double]$expectedRow.support) 0 "$($case.Name) support node $node"
            foreach ($suffix in @('_numeric_valid', '_display_valid', '_lambda_valid')) {
                Assert-Close $grid.Properties[$prefix + $suffix][$node] 1 0 "$($case.Name) $suffix node $node"
            }
            $maxRawError = [Math]::Max($maxRawError, [Math]::Abs($actualRaw - $expectedRaw))
            $maxVarianceError = [Math]::Max($maxVarianceError, [Math]::Abs($actualVariance - $expectedVariance))
            $maxLambdaError = [Math]::Max($maxLambdaError, [Math]::Abs($actualLambda - $expectedLambda))
            $maxEssError = [Math]::Max($maxEssError, [Math]::Abs($actualEss - $expectedEss))
            $exampleChecks.Add([PSCustomObject]@{ parameter_file = $case.Name; property = $prefix; weight_mode = $case.WeightMode; order = $case.Order; node_id = $node; expected_raw = $expectedRaw; actual_raw = $actualRaw; expected_kappa2 = $expectedVariance; actual_kappa2 = $actualVariance; expected_lambda = $expectedLambda; actual_lambda = $actualLambda; expected_ess = $expectedEss; actual_ess = $actualEss; support = [int]$expectedRow.support })
        }
        $checkedNodes = 9
        $referenceDescription = 'expected_window.csv: phi, matching weight_mode and order'
    }
    $exampleResults += [PSCustomObject]@{
        parameter_file = $case.Name; parameter_sha256 = $case.ParameterSha256; parameter_values_changed = $false;
        algorithm = $case.Algorithm; algorithm_selector_key = $case.CalculatorKey; calculator = $case.Calculator;
        weight_mode = $case.WeightMode; order = $case.Order; output_grid = $case.Grid; output_property = $prefix;
        export_file = [IO.Path]::GetFileName($case.ExportPath); verified_nodes = $checkedNodes; status = 'passed';
        reference = $referenceDescription; max_abs_raw_error = $maxRawError;
        max_abs_kappa2_error = $(if ($case.Algorithm -eq 'HOSC_WINDOW') { $maxVarianceError } else { $null });
        max_abs_lambda_error = $(if ($case.Algorithm -eq 'HOSC_WINDOW') { $maxLambdaError } else { $null });
        max_abs_ess_error = $(if ($case.Algorithm -eq 'HOSC_WINDOW') { $maxEssError } else { $null })
    }
}
$exampleChecks | Export-Csv -LiteralPath (Join-Path $EvidenceRoot 'gui_example_result_checks.csv') -NoTypeInformation -Encoding UTF8
$summary = [ordered]@{
    status = 'passed'; host_process_id = $hostPid; host_exit_code = $exitCode; isolated_host = $SGeMSRoot;
    template_invalid_diagnostic_fields = @('template_sampling', 'integer_direction_vectors');
    fixed_policy_rejections = @($policyCases | Select-Object Name,Algorithm,Field);
    template_valid_nodes = 27; template_support_min = 9; template_support_max = 75;
    window_invalid_diagnostic_field = 'windowl'; window_valid_rows = 9;
    same_process_recovery = $true; command_mode_no_modal_hang = $true;
    gui_highlight_visually_verified = $false;
    gui_parameter_files_executed = $exampleCases.Count; gui_parameter_values_unchanged = $true;
    gui_parameter_examples = $exampleResults;
    expected_window_sha256 = (Get-FileHash -LiteralPath $ExpectedWindowPath -Algorithm SHA256).Hash.ToLowerInvariant();
    gui_export_precision = 'IEEE-754 binary32';
    gui_export_tolerances = [ordered]@{ raw = 'max(1e-12, 2e-7*abs(reference))'; kappa2 = 'max(1e-10, 2e-7*abs(reference))'; lambda = 'max(1e-6, 2e-7*abs(reference))'; ess = 1e-5; support_and_masks = 0; template_raw = 1e-9 };
    host_sha256 = (Get-FileHash -LiteralPath $sgemsExe -Algorithm SHA256).Hash;
    template_dll_sha256 = (Get-FileHash -LiteralPath (Join-Path $SGeMSRoot 'plugins\Geostat\HOSC2008.dll') -Algorithm SHA256).Hash;
    window_dll_sha256 = (Get-FileHash -LiteralPath (Join-Path $SGeMSRoot 'plugins\Geostat\HOSCWindow.dll') -Algorithm SHA256).Hash
}
$summary | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $EvidenceRoot 'RESULT.json') -Encoding UTF8
Write-Host "PASS: recovery regression and all four unchanged GUI parameter examples, including two weighted sixth-order cases. Evidence: $EvidenceRoot"
