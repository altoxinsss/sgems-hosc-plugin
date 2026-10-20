#Requires -Version 5.1
[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$StageDir,
    [Parameter(Mandatory=$true)][string]$SGeMSRoot
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'host_common.ps1')
$StageDir = (Get-Item -LiteralPath $StageDir).FullName
$hostInfo = Resolve-HoscHost $SGeMSRoot
$prePath = Join-Path $StageDir 'source_manifest_prebuild.json'
$postPath = Join-Path $StageDir 'source_manifest_postbuild.json'
if ((Get-HoscSha256 $prePath) -ne (Get-HoscSha256 $postPath)) {
    throw 'The build inputs changed during compilation.'
}
$buildEntries = Get-Content -LiteralPath (Join-Path $StageDir 'build_manifest.json') -Raw | ConvertFrom-Json
$artifactNames = @('HOSC2008.dll','HOSCWindow.dll','HOSC_TEMPLATE.ui','HOSC_WINDOW.ui')
foreach ($name in $artifactNames) {
    $artifactEntries = @($buildEntries | Where-Object { [IO.Path]::GetFileName($_.path) -ceq $name })
    if ($artifactEntries.Count -ne 1 -or (Get-HoscSha256 (Join-Path $StageDir $name)) -ne $artifactEntries[0].sha256) {
        throw "Build artifact checksum mismatch: $name"
    }
}
$allSourceFiles = Get-Content -LiteralPath $prePath -Raw | ConvertFrom-Json
$sourceFiles = @($allSourceFiles | Where-Object {
    $_.relative_path -match '^(compat|shared|Template|MoveWindow|ui|tests)[\\/]' -and
    $_.relative_path -match '\.(cpp|h|ui)$'
})
foreach ($entry in $sourceFiles) {
    if ((Get-HoscSha256 (Join-Path $repo $entry.relative_path)) -ne $entry.sha256) {
        throw "Current source differs from the compiled source: $($entry.relative_path). Rebuild first."
    }
}
$targetDir = Join-Path $repo 'bin\win32'
New-Item -ItemType Directory -Path $targetDir -Force | Out-Null
foreach ($name in $artifactNames) { Copy-Item -LiteralPath (Join-Path $StageDir $name) -Destination $targetDir -Force }
$manifest = [ordered]@{
    schema_version = 1
    package_version = '2026.10.20'
    built_at_utc = [DateTime]::UtcNow.ToString('o')
    files = @($artifactNames | ForEach-Object { [ordered]@{name=$_;sha256=(Get-HoscSha256 (Join-Path $targetDir $_))} })
    runtime = [ordered]@{
        machine='x86'; compiler='VC7.1'; qt='3.3.5'
        required_files=@('MSVCR71.dll','MSVCP71.dll','qt-mt335.dll' | ForEach-Object { [ordered]@{name=$_} })
    }
    tested_host = @('sgems.exe','qt-mt335.dll','GsTLAppli_grid.dll','GsTLAppli_utils.dll','MSVCR71.dll','MSVCP71.dll' | ForEach-Object {
        [ordered]@{name=$_;sha256=(Get-HoscSha256 (Join-Path $hostInfo.Root $_))}
    })
    source_manifest = 'build_source_manifest.json'
    source_manifest_note = 'Compiled C++/header/UI inputs, verified against this checkout; host hashes identify the tested distribution.'
}
Write-HoscJson (Join-Path $targetDir 'build_source_manifest.json') @($sourceFiles | ForEach-Object {
    [ordered]@{path=$_.relative_path.Replace('\','/');bytes=$_.bytes;sha256=$_.sha256}
})
Write-HoscJson (Join-Path $targetDir 'manifest.json') $manifest
$package = Get-HoscPackage $targetDir
Assert-HoscCompatibility $hostInfo $package
Write-Host "Release artifacts and source hashes verified: $targetDir"
Write-Host 'Run Install.cmd and RunExamples.cmd against a separate compatible host before publishing.'
