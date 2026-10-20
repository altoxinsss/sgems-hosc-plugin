#Requires -Version 5.1
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$BackupManifest,
    [string]$SGeMSRoot = $env:HOSC_SGEMS_RUNTIME,
    [switch]$CheckOnly
)
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'host_common.ps1')
$installLock = $null
$changed = New-Object 'System.Collections.Generic.List[string]'
$snapshotDir = $null
try {
    $hostInfo = Resolve-HoscHost $SGeMSRoot
    Assert-HoscHostStopped $hostInfo
    $manifestPath = Get-HoscFullPath $BackupManifest
    $backupDir = Split-Path -Parent $manifestPath
    if ((Split-Path -Leaf $manifestPath) -ne 'backup_manifest.json' -or (Get-HoscFullPath (Split-Path -Parent $backupDir)) -ne $hostInfo.BackupRoot -or (Split-Path -Leaf $backupDir) -notmatch '^hosc_[0-9]{8}_[0-9]{6}_[0-9a-f]{32}$') {
        throw 'Select an explicit backup_manifest.json created by this installer under this host plugin_backups folder.'
    }
    Assert-HoscPlainDirectory $hostInfo.BackupRoot
    Assert-HoscPlainDirectory $backupDir
    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    if ($manifest.schema_version -ne 1 -or $manifest.kind -ne 'hosc-plugin-backup' -or (Get-HoscFullPath $manifest.host_root) -ne $hostInfo.Root) { throw 'Backup manifest does not belong to this SGeMS installation.' }
    $records = @($manifest.files)
    if ($records.Count -ne 4) { throw 'Invalid backup file count.' }
    $alreadyRestored = $true
    foreach ($name in $script:HoscArtifactNames) {
        $matches = @($records | Where-Object { $_.name -ceq $name })
        if ($matches.Count -ne 1) { throw "Invalid backup entry: $name" }
        $entry = $matches[0]
        if ($entry.existed -isnot [bool] -or $entry.installed_sha256 -notmatch '^[0-9a-fA-F]{64}$') { throw "Invalid backup metadata: $name" }
        if ($entry.existed) {
            $source = Join-Path $backupDir $name
            Assert-HoscTargetFile $source
            if ($entry.original_sha256 -notmatch '^[0-9a-fA-F]{64}$' -or -not (Test-Path -LiteralPath $source -PathType Leaf) -or (Get-HoscSha256 $source) -ne $entry.original_sha256) { throw "Original backup checksum mismatch: $name" }
        }
        $target = Join-Path $hostInfo.PluginDir $name
        Assert-HoscTargetFile $target
        $exists = Test-Path -LiteralPath $target -PathType Leaf
        $hash = $null
        if ($exists) { $hash = Get-HoscSha256 $target }
        $matchesOriginal = (($entry.existed -and $exists -and $hash -eq $entry.original_sha256) -or (-not $entry.existed -and -not $exists))
        if (-not $matchesOriginal) { $alreadyRestored = $false }
        if (-not $matchesOriginal -and (-not $exists -or $hash -ne $entry.installed_sha256)) { throw "Refusing to overwrite a file changed after installation: $name" }
    }
    if ($CheckOnly) { Write-Host "PASS: backup integrity and restoration targets. Already restored: $alreadyRestored"; exit 0 }
    $installLock = Enter-HoscInstallLock $hostInfo
    if ($alreadyRestored) { Write-Host 'The previous file state is already restored. No files were changed.'; exit 0 }
    $snapshotDir = New-HoscBackupDirectory $hostInfo 'uninstall'
    $snapshot = @()
    foreach ($entry in $records) {
        $target = Join-Path $hostInfo.PluginDir $entry.name
        $exists = Test-Path -LiteralPath $target -PathType Leaf
        $hash = $null
        if ($exists) {
            $hash = Get-HoscSha256 $target
            [IO.File]::Copy($target, (Join-Path $snapshotDir $entry.name), $false)
            if ((Get-HoscSha256 (Join-Path $snapshotDir $entry.name)) -ne $hash) { throw 'Restoration recovery-copy verification failed.' }
        }
        $snapshot += [PSCustomObject]@{ name = $entry.name; existed = [bool]$exists; sha256 = $hash }
    }
    Write-HoscJson (Join-Path $snapshotDir 'restoration_snapshot.json') ([PSCustomObject]@{ host_root = $hostInfo.Root; backup_manifest = $manifestPath; files = $snapshot })
    Assert-HoscHostStopped $hostInfo
    foreach ($entry in $records) {
        $target = Join-Path $hostInfo.PluginDir $entry.name
        Assert-HoscTargetFile $target
        $before = @($snapshot | Where-Object { $_.name -eq $entry.name })[0]
        if ($before.existed) {
            if (-not (Test-Path -LiteralPath $target -PathType Leaf) -or (Get-HoscSha256 $target) -ne $before.sha256) { throw "Target changed during restoration: $($entry.name)" }
        }
        elseif (Test-Path -LiteralPath $target) { throw "Target appeared during restoration: $($entry.name)" }
        if ($entry.existed) { Set-HoscFileAtomic (Join-Path $backupDir $entry.name) $target }
        elseif (Test-Path -LiteralPath $target -PathType Leaf) { [IO.File]::Delete($target) }
        $changed.Add($entry.name)
    }
    $manifest.status = 'restored'
    Write-HoscJson $manifestPath $manifest
    Write-Host 'Previous DLL/UI file state restored. Other plugins were not changed. Backups remain available.'
}
catch {
    $failureMessage = $_.Exception.Message
    if ($changed.Count -gt 0) {
        Write-Warning 'Restoration failed. Returning changed files to their state before this operation.'
        $reverseChanged = @($changed.ToArray())
        [Array]::Reverse($reverseChanged)
        foreach ($name in $reverseChanged) {
            try {
                $record = @($snapshot | Where-Object { $_.name -eq $name })[0]
                $original = @($records | Where-Object { $_.name -eq $name })[0]
                $target = Join-Path $hostInfo.PluginDir $name
                Assert-HoscTargetFile $target
                if ($original.existed) {
                    if (-not (Test-Path -LiteralPath $target -PathType Leaf) -or (Get-HoscSha256 $target) -ne $original.original_sha256) { throw 'Target changed externally; preserving it.' }
                }
                elseif (Test-Path -LiteralPath $target) { throw 'Target appeared externally; preserving it.' }
                if ($record.existed) { Set-HoscFileAtomic (Join-Path $snapshotDir $name) $target }
                elseif (Test-Path -LiteralPath $target -PathType Leaf) { [IO.File]::Delete($target) }
            }
            catch { Write-Warning "Could not recover $name : $($_.Exception.Message)" }
        }
    }
    Write-Error -Message $failureMessage -ErrorAction Continue
    if ($snapshotDir) { Write-Host "Recovery snapshot retained: $snapshotDir" }
    exit 1
}
finally { if ($installLock) { $installLock.ReleaseMutex(); $installLock.Dispose() } }
