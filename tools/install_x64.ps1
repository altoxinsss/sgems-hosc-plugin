param(
    [Parameter(Mandatory=$true)][string]$TargetRoot,
    [switch]$VerifyOnly
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$packageRoot = Join-Path $repoRoot 'bin\win64'
$hostRoot = (Resolve-Path -LiteralPath $TargetRoot).Path
$hostExe = Join-Path $hostRoot 'sgems-x64.exe'
if (!(Test-Path -LiteralPath $hostExe -PathType Leaf)) { throw "Missing sgems-x64.exe: $hostRoot" }
$bytes = [IO.File]::ReadAllBytes($hostExe)
$peOffset = [BitConverter]::ToInt32($bytes, 0x3c)
if ([BitConverter]::ToUInt32($bytes, $peOffset) -ne 0x4550 -or
    [BitConverter]::ToUInt16($bytes, $peOffset + 4) -ne 0x8664) { throw 'Target is not an x64 PE host' }
$hostRows = @(Import-Csv -LiteralPath (Join-Path $packageRoot 'host_libraries.csv'))
$packageRows = @(Import-Csv -LiteralPath (Join-Path $packageRoot 'SHA256SUMS.csv'))
$expected = @('HOSC2008.dll', 'HOSCWindow.dll', 'HOSC_TEMPLATE.ui', 'HOSC_WINDOW.ui')
if ($packageRows.Count -ne 4 -or @(Compare-Object ($packageRows.file | Sort-Object) ($expected | Sort-Object)).Count) {
    throw 'Unexpected release manifest'
}
foreach ($row in $hostRows) {
    if ([IO.Path]::GetFileName($row.file) -ne $row.file) { throw 'Invalid host manifest path' }
    $hostFile = Join-Path $hostRoot $row.file
    if (!(Test-Path -LiteralPath $hostFile -PathType Leaf) -or
        (Get-FileHash -LiteralPath $hostFile -Algorithm SHA256).Hash -ne $row.sha256) {
        throw "Unsupported host library: $($row.file). No files were installed."
    }
}
foreach ($row in $packageRows) {
    $releaseFile = Join-Path $packageRoot $row.file
    if (!(Test-Path -LiteralPath $releaseFile -PathType Leaf) -or
        (Get-FileHash -LiteralPath $releaseFile -Algorithm SHA256).Hash -ne $row.sha256) {
        throw "Release checksum mismatch: $($row.file)"
    }
}
Write-Output "Verified x64 host, $($hostRows.Count) host libraries and 4 release files."
if ($VerifyOnly) { return }
$running = @(Get-CimInstance Win32_Process -Filter "Name='sgems-x64.exe'" |
    Where-Object { !$_.ExecutablePath -or $_.ExecutablePath -ieq $hostExe })
if ($running.Count) { throw 'Close this SGeMS host before installation.' }
$destination = Join-Path $hostRoot 'plugins\Geostat'
if (!(Test-Path -LiteralPath $destination -PathType Container)) { throw "Missing plugin directory: $destination" }
$backup = Join-Path $hostRoot ('plugin_backups\hosc_x64_' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff'))
New-Item -ItemType Directory -Path $backup | Out-Null
$records = @()
foreach ($row in $packageRows) {
    $target = Join-Path $destination $row.file
    $exists = Test-Path -LiteralPath $target -PathType Leaf
    $oldHash = ''
    if ($exists) {
        $oldHash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
        Copy-Item -LiteralPath $target -Destination (Join-Path $backup $row.file)
        if ((Get-FileHash -LiteralPath (Join-Path $backup $row.file) -Algorithm SHA256).Hash -ne $oldHash) {
            throw "Backup verification failed: $target"
        }
    }
    $records += [pscustomobject]@{file=$row.file;previously_present=$exists;previous_sha256=$oldHash;new_sha256=$row.sha256}
}
$records | Export-Csv -LiteralPath (Join-Path $backup 'installation.csv') -NoTypeInformation -Encoding UTF8
try {
    foreach ($row in $packageRows) {
        $target = Join-Path $destination $row.file
        Copy-Item -LiteralPath (Join-Path $packageRoot $row.file) -Destination $target -Force
        if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -ne $row.sha256) { throw "Installed checksum mismatch: $target" }
    }
} catch {
    Write-Warning "Installation incomplete. Verified previous files are in $backup; no host libraries were changed."
    throw
}
Write-Output "Installed four plugin files in $destination"
Write-Output "Previous plugin files and installation record: $backup"
