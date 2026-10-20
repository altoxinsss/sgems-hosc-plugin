#Requires -Version 5.1
[CmdletBinding()]
param(
    [Parameter(Position=0)][string]$SGeMSRoot = $env:HOSC_SGEMS_RUNTIME,
    [string]$OutputDirectory = '',
    [string]$WorkDirectory = '',
    [int]$TimeoutSeconds = 120
)
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'host_common.ps1')
$repo = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($SGeMSRoot)) { $SGeMSRoot = Read-Host 'Enter the SGeMS folder containing sgems.exe' }
$hostInfo = Resolve-HoscHost $SGeMSRoot
$package = Get-HoscPackage (Join-Path $repo 'bin\win32')
Assert-HoscCompatibility $hostInfo $package
foreach ($entry in $package.Files) {
    $installed = Join-Path $hostInfo.PluginDir $entry.name
    if (!(Test-Path -LiteralPath $installed -PathType Leaf) -or (Get-HoscSha256 $installed) -ne $entry.sha256) {
        throw "Install the current package first (Install.cmd). Installed file differs: $($entry.name)"
    }
}
if (!$OutputDirectory) { $OutputDirectory = Join-Path $repo 'reviewer-output' }
$OutputDirectory = [IO.Path]::GetFullPath($OutputDirectory)
$runName = 'run_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [Guid]::NewGuid().ToString('N').Substring(0,8)
$result = Join-Path $OutputDirectory $runName
if (!$WorkDirectory) { $WorkDirectory = [IO.Path]::GetTempPath() }
if (!(Test-Path -LiteralPath $WorkDirectory -PathType Container)) { throw 'WorkDirectory must be an existing writable folder.' }
$WorkDirectory = [IO.Path]::GetFullPath($WorkDirectory)
# SGeMS 2.0 command-file loading uses the legacy ANSI API. Repository and
# output paths may contain Unicode; the independent run uses an ASCII path.
if ($WorkDirectory -match '[^\x00-\x7F]') {
    try {
        $filesystem = New-Object -ComObject Scripting.FileSystemObject
        $WorkDirectory = $filesystem.GetFolder($WorkDirectory).ShortPath
    }
    catch { }
}
if ($WorkDirectory -match '[^\x00-\x7F]') {
    throw 'Choose an existing ASCII-only writable folder with -WorkDirectory, for example C:\HOSC-work.'
}
$work = Join-Path $WorkDirectory ('hosc_' + [Guid]::NewGuid().ToString('N'))
$runtime = Join-Path $work 'runtime'
$evidence = Join-Path $work 'results'
New-Item -ItemType Directory -Path $runtime -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $runtime 'plugins\Geostat') -Force | Out-Null
Get-ChildItem -LiteralPath $hostInfo.Root -Filter '*.dll' -File | Copy-Item -Destination $runtime
Copy-Item -LiteralPath $hostInfo.Exe -Destination (Join-Path $runtime 'sgems.exe')
# Some SGeMS distributions include their Python 2.4 standard library here.
# Preserve it in the isolated host when available; no external Python is used.
foreach ($folderName in @('lib', 'encodings')) {
    $folder = Join-Path $hostInfo.Root $folderName
    if (Test-Path -LiteralPath $folder -PathType Container) {
        Assert-HoscPlainDirectory $folder
        Copy-Item -LiteralPath $folder -Destination $runtime -Recurse
    }
}
foreach ($entry in $package.Files) {
    Copy-Item -LiteralPath (Join-Path $hostInfo.PluginDir $entry.name) -Destination (Join-Path $runtime ('plugins\Geostat\' + $entry.name))
}
$fixture = Join-Path $work 'selftest_3d_cartesian.sgems'
$reference = Join-Path $work 'expected_window.csv'
$examples = Join-Path $work 'examples'
Copy-Item -LiteralPath (Join-Path $repo 'tests\fixtures\manual\selftest_3d_cartesian.sgems') -Destination $fixture
Copy-Item -LiteralPath (Join-Path $repo 'tests\fixtures\manual\expected_window.csv') -Destination $reference
New-Item -ItemType Directory -Path $examples | Out-Null
foreach ($name in @('template_3d.par', 'window_uniform.par', 'window_coordinate_gaussian.par', 'window_soft_high_porosity.par')) {
    Copy-Item -LiteralPath (Join-Path $repo ('examples\' + $name)) -Destination (Join-Path $examples $name)
}
try {
    & (Join-Path $PSScriptRoot 'verify_runtime_examples.ps1') -SGeMSRoot $runtime -FixturePath $fixture -ExpectedWindowPath $reference -ExampleDirectory $examples -EvidenceRoot $evidence -TimeoutSeconds $TimeoutSeconds
    New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
    Copy-Item -LiteralPath $evidence -Destination $result -Recurse
    Copy-Item -LiteralPath $fixture -Destination (Join-Path $result 'input.sgems')
    Write-Host 'PASS: four unchanged GUI PAR files (template, uniform r4, Gaussian r6, soft high-porosity r6), values, masks, support, ESS, and invalid-input recovery.'
    Write-Host "Results: $result"
    Write-Host 'Load input.sgems, template_3d.sgems, window_uniform.sgems, window_coordinate_gaussian.sgems, and window_soft_high_porosity.sgems to inspect the examples.'
    Write-Host "Independent runtime and diagnostic files: $work"
}
catch {
    Write-Host "Example run did not pass. Diagnostic files are retained in: $work"
    throw
}
