#Requires -Version 5.1
[CmdletBinding()]
param(
    [Parameter(Position = 0)][string]$SGeMSRoot = $env:HOSC_SGEMS_RUNTIME,
    [string]$PackageDir = '',
    [switch]$CheckOnly
)
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'host_common.ps1')
$installLock = $null
$backupDir = $null
$backupManifest = $null
$changed = New-Object 'System.Collections.Generic.List[string]'
try {
    if ([string]::IsNullOrWhiteSpace($SGeMSRoot)) {
        $SGeMSRoot = Read-Host 'Enter the SGeMS 2.0 folder containing sgems.exe'
    }
    if ([string]::IsNullOrWhiteSpace($PackageDir)) { $PackageDir = Join-Path (Split-Path -Parent $PSScriptRoot) 'bin\win32' }
    $hostInfo = Resolve-HoscHost $SGeMSRoot
    $package = Get-HoscPackage $PackageDir
    Assert-HoscCompatibility $hostInfo $package
    Assert-HoscHostStopped $hostInfo
    $allInstalled = $true
    foreach ($entry in $package.Files) {
        $target = Join-Path $hostInfo.PluginDir $entry.name
        Assert-HoscTargetFile $target
        if (-not (Test-Path -LiteralPath $target -PathType Leaf) -or (Get-HoscSha256 $target) -ne $entry.sha256) { $allInstalled = $false }
    }
    if ($CheckOnly) {
        Write-Host 'PASS: package checksums, x86 architecture, VC7.1/Qt runtime, target paths, and stopped host.'
        Write-Host "SGeMS folder: $($hostInfo.Root)"
        Write-Host "Packaged plugin already installed: $allInstalled"
        exit 0
    }
    $installLock = Enter-HoscInstallLock $hostInfo
    if ($allInstalled) {
        Write-Host 'The exact packaged DLL/UI files are already installed. No files were changed.'
        exit 0
    }
    $backupDir = New-HoscBackupDirectory $hostInfo 'hosc'
    $records = @()
    foreach ($entry in $package.Files) {
        $target = Join-Path $hostInfo.PluginDir $entry.name
        $existed = Test-Path -LiteralPath $target -PathType Leaf
        $oldHash = $null
        if ($existed) {
            $oldHash = Get-HoscSha256 $target
            [IO.File]::Copy($target, (Join-Path $backupDir $entry.name), $false)
            if ((Get-HoscSha256 (Join-Path $backupDir $entry.name)) -ne $oldHash) { throw "Backup verification failed: $($entry.name)" }
        }
        $records += [PSCustomObject]@{ name = $entry.name; existed = [bool]$existed; original_sha256 = $oldHash; installed_sha256 = $entry.sha256.ToLowerInvariant() }
    }
    $backupManifest = [PSCustomObject]@{
        schema_version = 1; kind = 'hosc-plugin-backup'; host_root = $hostInfo.Root
        created_utc = [DateTime]::UtcNow.ToString('o'); status = 'prepared'
        package_manifest_sha256 = Get-HoscSha256 $package.ManifestPath; files = $records
    }
    $backupManifestPath = Join-Path $backupDir 'backup_manifest.json'
    Write-HoscJson $backupManifestPath $backupManifest
    Assert-HoscHostStopped $hostInfo
    foreach ($entry in $package.Files) {
        $target = Join-Path $hostInfo.PluginDir $entry.name
        $record = @($records | Where-Object { $_.name -eq $entry.name })[0]
        Assert-HoscTargetFile $target
        if ($record.existed) {
            if (-not (Test-Path -LiteralPath $target -PathType Leaf) -or (Get-HoscSha256 $target) -ne $record.original_sha256) { throw "Target changed during installation: $($entry.name)" }
        }
        elseif (Test-Path -LiteralPath $target) { throw "Target appeared during installation: $($entry.name)" }
        $source = Join-Path $package.Directory $entry.name
        if ((Get-HoscSha256 $source) -ne $entry.sha256) { throw "Package changed during installation: $($entry.name)" }
        Set-HoscFileAtomic $source $target
        $changed.Add($entry.name)
        if ((Get-HoscSha256 $target) -ne $entry.sha256) { throw "Installed checksum mismatch: $($entry.name)" }
    }
    $backupManifest.status = 'installed'
    Write-HoscJson $backupManifestPath $backupManifest
    Write-Host 'Installation complete. Start SGeMS and select HOSC_TEMPLATE or HOSC_WINDOW.'
    Write-Host "Backup manifest: $backupManifestPath"
    Write-Host 'To restore the previous files, close SGeMS and run:'
    Write-Host ('powershell.exe -NoProfile -ExecutionPolicy Bypass -File "' + (Join-Path $PSScriptRoot 'uninstall_plugin.ps1') + '" -SGeMSRoot "' + $hostInfo.Root + '" -BackupManifest "' + $backupManifestPath + '"')
}
catch {
    $failureMessage = $_.Exception.Message
    $rollbackFailed = $false
    if ($changed.Count -gt 0 -and $backupManifest) {
        Write-Warning 'Installation failed. Restoring files changed by this installation.'
        $reverseChanged = @($changed.ToArray())
        [Array]::Reverse($reverseChanged)
        foreach ($name in $reverseChanged) {
            try {
                $record = @($backupManifest.files | Where-Object { $_.name -eq $name })[0]
                $target = Join-Path $hostInfo.PluginDir $name
                Assert-HoscTargetFile $target
                if (-not (Test-Path -LiteralPath $target -PathType Leaf) -or (Get-HoscSha256 $target) -ne $record.installed_sha256) { throw 'Target changed externally; preserving it.' }
                if ($record.existed) { Set-HoscFileAtomic (Join-Path $backupDir $name) $target }
                else { [IO.File]::Delete($target) }
            }
            catch { $rollbackFailed = $true; Write-Warning "Could not roll back $name : $($_.Exception.Message)" }
        }
        if ($rollbackFailed) { $backupManifest.status = 'rollback_incomplete' } else { $backupManifest.status = 'rolled_back' }
        try { Write-HoscJson (Join-Path $backupDir 'backup_manifest.json') $backupManifest } catch { Write-Warning 'Could not update backup status.' }
    }
    Write-Error -Message $failureMessage -ErrorAction Continue
    if ($backupDir) { Write-Host "Recovery backup retained: $backupDir" }
    exit 1
}
finally {
    if ($installLock) { $installLock.ReleaseMutex(); $installLock.Dispose() }
}
