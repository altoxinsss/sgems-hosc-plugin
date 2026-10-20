#Requires -Version 5.1
# Shared, compiler-free checks for the packaged 32-bit SGeMS plugins.
$script:HoscArtifactNames = @('HOSC2008.dll', 'HOSCWindow.dll', 'HOSC_TEMPLATE.ui', 'HOSC_WINDOW.ui')

function Get-HoscSha256([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Get-HoscFullPath([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if ($full -eq [IO.Path]::GetPathRoot($full)) { return $full }
    return $full.TrimEnd([IO.Path]::DirectorySeparatorChar)
}

function Assert-HoscPlainDirectory([string]$Path) {
    $item = Get-Item -LiteralPath $Path -ErrorAction Stop
    if (-not $item.PSIsContainer -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
        throw "Expected an ordinary directory (not a junction or symbolic link): $Path"
    }
}

function Resolve-HoscHost([string]$SGeMSRoot) {
    if ([string]::IsNullOrWhiteSpace($SGeMSRoot)) {
        throw 'Specify the folder containing sgems.exe with -SGeMSRoot.'
    }
    $root = Get-HoscFullPath $SGeMSRoot.Trim().Trim('"')
    Assert-HoscPlainDirectory $root
    $exe = Join-Path $root 'sgems.exe'
    if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
        throw "sgems.exe was not found in: $root"
    }
    $plugins = Join-Path $root 'plugins'
    $geostat = Join-Path $plugins 'Geostat'
    Assert-HoscPlainDirectory $plugins
    Assert-HoscPlainDirectory $geostat
    $backupRoot = Join-Path $root 'plugin_backups'
    if (Test-Path -LiteralPath $backupRoot) { Assert-HoscPlainDirectory $backupRoot }
    return [PSCustomObject]@{ Root = $root; Exe = $exe; PluginDir = $geostat; BackupRoot = $backupRoot }
}

function Assert-HoscHostStopped($HostInfo) {
    foreach ($process in @(Get-Process -Name 'sgems' -ErrorAction SilentlyContinue)) {
        $processPath = $null
        try { $processPath = $process.Path } catch { }
        if ([string]::IsNullOrWhiteSpace($processPath)) {
            throw "A running sgems.exe process (PID $($process.Id)) cannot be identified. Close SGeMS before installation or restoration."
        }
        if ((Get-HoscFullPath $processPath) -eq $HostInfo.Exe) {
            throw "SGeMS is running from the selected folder (PID $($process.Id)). Save the project and close this SGeMS instance first."
        }
    }
}

function Get-HoscPeInfo([string]$Path) {
    $stream = [IO.File]::Open($Path, [IO.FileMode]::Open, [IO.FileAccess]::Read, [IO.FileShare]::Read)
    $reader = New-Object IO.BinaryReader($stream)
    try {
        if ($stream.Length -lt 64 -or $reader.ReadUInt16() -ne 0x5a4d) { throw 'Missing DOS header.' }
        $stream.Position = 0x3c
        $peOffset = $reader.ReadUInt32()
        if ($peOffset + 24 -gt $stream.Length) { throw 'Invalid PE header offset.' }
        $stream.Position = $peOffset
        if ($reader.ReadUInt32() -ne 0x4550) { throw 'Missing PE signature.' }
        $machine = $reader.ReadUInt16()
        $sectionCount = $reader.ReadUInt16()
        $stream.Position = $peOffset + 20
        $optionalSize = $reader.ReadUInt16()
        $optionalOffset = $peOffset + 24
        if ($optionalOffset + $optionalSize + 40 * $sectionCount -gt $stream.Length) { throw 'Truncated PE headers.' }
        $stream.Position = $optionalOffset
        $magic = $reader.ReadUInt16()
        if ($magic -eq 0x10b) { $directoryOffset = 96 }
        elseif ($magic -eq 0x20b) { $directoryOffset = 112 }
        else { throw 'Unknown PE optional header.' }
        if ($optionalSize -lt $directoryOffset + 16) { throw 'Missing PE import directory.' }
        $stream.Position = $optionalOffset + $directoryOffset + 8
        $importRva = $reader.ReadUInt32()
        $importSize = $reader.ReadUInt32()
        $sections = @()
        for ($index = 0; $index -lt $sectionCount; $index++) {
            $stream.Position = $optionalOffset + $optionalSize + 40 * $index + 8
            $virtualSize = $reader.ReadUInt32()
            $virtualAddress = $reader.ReadUInt32()
            $rawSize = $reader.ReadUInt32()
            $rawOffset = $reader.ReadUInt32()
            $sections += [PSCustomObject]@{ Rva = $virtualAddress; Size = [Math]::Max($virtualSize, $rawSize); RawSize = $rawSize; Raw = $rawOffset }
        }
        $mapRva = {
            param([uint32]$Rva)
            foreach ($section in $sections) {
                if ([int64]$Rva -ge $section.Rva -and [int64]$Rva -lt ([int64]$section.Rva + $section.Size)) {
                    $delta = [int64]$Rva - $section.Rva
                    if ($delta -ge $section.RawSize -or $section.Raw + $delta -ge $stream.Length) { throw 'PE RVA has no file data.' }
                    return [int64]($section.Raw + $delta)
                }
            }
            if ($Rva -lt $optionalOffset + $optionalSize) { return [int64]$Rva }
            throw 'PE RVA is outside the image sections.'
        }
        $imports = @()
        if ($importRva -ne 0) {
            $importOffset = & $mapRva $importRva
            $terminated = $false
            for ($index = 0; $index -lt 2048; $index++) {
                $entryOffset = $importOffset + 20 * $index
                if ($entryOffset + 20 -gt $stream.Length) { throw 'Truncated PE import table.' }
                $stream.Position = $entryOffset
                $entry = @($reader.ReadUInt32(), $reader.ReadUInt32(), $reader.ReadUInt32(), $reader.ReadUInt32(), $reader.ReadUInt32())
                if (($entry | Where-Object { $_ -ne 0 }).Count -eq 0) { $terminated = $true; break }
                $stream.Position = & $mapRva $entry[3]
                $bytes = New-Object 'System.Collections.Generic.List[byte]'
                $nameEnded = $false
                for ($letter = 0; $letter -lt 512 -and $stream.Position -lt $stream.Length; $letter++) {
                    $value = $reader.ReadByte()
                    if ($value -eq 0) { $nameEnded = $true; break }
                    $bytes.Add($value)
                }
                if (-not $nameEnded) { throw 'Invalid PE import name.' }
                $imports += [Text.Encoding]::ASCII.GetString($bytes.ToArray())
            }
            if (-not $terminated) { throw 'PE import table has no terminator.' }
        }
        return [PSCustomObject]@{ Machine = $machine; Imports = @($imports | Sort-Object -Unique) }
    }
    catch { throw "Cannot inspect PE file '$Path': $($_.Exception.Message)" }
    finally { $reader.Dispose(); $stream.Dispose() }
}

function Assert-HoscLeafName([string]$Name) {
    if ([string]::IsNullOrWhiteSpace($Name) -or $Name -ne [IO.Path]::GetFileName($Name) -or $Name -match '[:/\\]' -or $Name -eq '.' -or $Name -eq '..') {
        throw "Invalid filename in manifest: $Name"
    }
}

function Get-HoscPackage([string]$PackageDir) {
    $directory = Get-HoscFullPath $PackageDir
    Assert-HoscPlainDirectory $directory
    $manifestPath = Join-Path $directory 'manifest.json'
    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    if ($manifest.schema_version -ne 1) { throw 'Unsupported package manifest schema_version (expected 1).' }
    if ($manifest.runtime.machine -ne 'x86') { throw 'The package manifest must identify the x86 runtime.' }
    $files = @($manifest.files)
    if ($files.Count -ne $script:HoscArtifactNames.Count) { throw 'Package manifest must list exactly the four HOSC DLL/UI artifacts.' }
    foreach ($name in $script:HoscArtifactNames) {
        $entry = @($files | Where-Object { $_.name -ceq $name })
        if ($entry.Count -ne 1 -or $entry[0].sha256 -notmatch '^[0-9a-fA-F]{64}$') { throw "Missing or invalid manifest entry: $name" }
        $path = Join-Path $directory $name
        if (-not (Test-Path -LiteralPath $path -PathType Leaf) -or (Get-HoscSha256 $path) -ne $entry[0].sha256) {
            throw "Package checksum mismatch: $name. Download a complete, unchanged release."
        }
    }
    return [PSCustomObject]@{ Directory = $directory; Manifest = $manifest; ManifestPath = $manifestPath; Files = $files }
}

function Assert-HoscCompatibility($HostInfo, $Package) {
    $hostPe = Get-HoscPeInfo $HostInfo.Exe
    if ($hostPe.Machine -ne 0x14c) { throw 'This package requires the 32-bit (x86) SGeMS 2.0 host.' }
    $runtimeNames = @('MSVCR71.dll', 'MSVCP71.dll', 'qt-mt335.dll')
    $systemImports = '^(KERNEL32|USER32|GDI32|ADVAPI32|SHELL32|OLE32|OLEAUT32|COMDLG32|COMCTL32|WS2_32|WSOCK32|WINMM|VERSION|RPCRT4|SHLWAPI|IMM32|NTDLL|OPENGL32|GLU32|WINSPOOL|MPR|NETAPI32|IPHLPAPI|WININET|CRYPT32|SECUR32|SETUPAPI|MSIMG32)\.dll$'
    $runtimeNames += @($hostPe.Imports | Where-Object { $_ -notmatch $systemImports })
    foreach ($name in @('HOSC2008.dll', 'HOSCWindow.dll')) {
        $pe = Get-HoscPeInfo (Join-Path $Package.Directory $name)
        if ($pe.Machine -ne 0x14c) { throw "Plugin is not x86: $name" }
        if ($pe.Imports -notcontains 'MSVCR71.dll' -or $pe.Imports -notcontains 'MSVCP71.dll' -or @($pe.Imports | Where-Object { $_ -match '^MSVC[RP](80|90|100|110|120|140)' }).Count -gt 0) {
            throw "Plugin does not match the VC7.1 runtime contract: $name"
        }
        foreach ($dependency in $pe.Imports) {
            Assert-HoscLeafName $dependency
            if ($dependency -match '^(GsTLAppli_|qt-|MSVC)') { $runtimeNames += $dependency }
            elseif ($dependency -notmatch $systemImports) {
                throw "Unexpected plugin dependency '$dependency' in $name; the release package needs review."
            }
        }
    }
    foreach ($entry in @($Package.Manifest.runtime.required_files)) {
        Assert-HoscLeafName $entry.name
        $runtimeNames += $entry.name
        $path = Join-Path $HostInfo.Root $entry.name
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Required SGeMS runtime file is missing: $path" }
        if ($entry.sha256) {
            if ($entry.sha256 -notmatch '^[0-9a-fA-F]{64}$' -or (Get-HoscSha256 $path) -ne $entry.sha256) {
                throw "SGeMS runtime differs from the tested release: $($entry.name). Use the supported SGeMS 2.0 distribution documented in README.md."
            }
        }
    }
    foreach ($name in @($runtimeNames | Sort-Object -Unique)) {
        Assert-HoscLeafName $name
        $path = Join-Path $HostInfo.Root $name
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Required SGeMS runtime file is missing: $path" }
        if ((Get-HoscPeInfo $path).Machine -ne 0x14c) { throw "Runtime file is not x86: $name" }
    }
    if (@($hostPe.Imports | Where-Object { $_ -match '^MSVC[RP](80|90|100|110|120|140)' }).Count -gt 0) {
        throw 'This SGeMS executable uses a different C++ runtime. The release supports the VC7.1/Qt 3.3.5 host only.'
    }
    $legacy = Join-Path $HostInfo.PluginDir 'HOSC_AVL.dll'
    if (Test-Path -LiteralPath $legacy) {
        Write-Warning 'An older HOSC_AVL.dll is present. It is not changed by this installer; review old HOSC registrations if duplicate algorithm names appear.'
    }
}

function Enter-HoscInstallLock($HostInfo) {
    $sha = [Security.Cryptography.SHA256]::Create()
    try { $key = [BitConverter]::ToString($sha.ComputeHash([Text.Encoding]::UTF8.GetBytes($HostInfo.Root.ToLowerInvariant()))).Replace('-', '') }
    finally { $sha.Dispose() }
    $mutex = New-Object Threading.Mutex($false, ('Local\HOSCInstall_' + $key))
    try { $acquired = $mutex.WaitOne(0) }
    catch [Threading.AbandonedMutexException] { $acquired = $true }
    if (-not $acquired) { $mutex.Dispose(); throw 'Another HOSC install or restoration is using this SGeMS folder.' }
    return $mutex
}

function New-HoscBackupDirectory($HostInfo, [string]$Prefix) {
    if (-not (Test-Path -LiteralPath $HostInfo.BackupRoot)) { [void][IO.Directory]::CreateDirectory($HostInfo.BackupRoot) }
    Assert-HoscPlainDirectory $HostInfo.BackupRoot
    $name = $Prefix + '_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [Guid]::NewGuid().ToString('N')
    $directory = Join-Path $HostInfo.BackupRoot $name
    [void][IO.Directory]::CreateDirectory($directory)
    return $directory
}

function Write-HoscJson([string]$Path, $Value) {
    $json = $Value | ConvertTo-Json -Depth 12
    [IO.File]::WriteAllText($Path, $json + [Environment]::NewLine, (New-Object Text.UTF8Encoding($false)))
}

function Set-HoscFileAtomic([string]$Source, [string]$Destination) {
    $pending = $Destination + '.' + [Guid]::NewGuid().ToString('N') + '.pending'
    try {
        [IO.File]::Copy($Source, $pending, $false)
        if ((Get-HoscSha256 $Source) -ne (Get-HoscSha256 $pending)) { throw "Copy verification failed: $Destination" }
        if ([IO.File]::Exists($Destination)) { [IO.File]::Replace($pending, $Destination, [System.Management.Automation.Language.NullString]::Value) }
        else { [IO.File]::Move($pending, $Destination) }
    }
    finally { if ([IO.File]::Exists($pending)) { [IO.File]::Delete($pending) } }
}

function Assert-HoscTargetFile([string]$Path) {
    if (Test-Path -LiteralPath $Path) {
        $item = Get-Item -LiteralPath $Path
        if ($item.PSIsContainer -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
            throw "Plugin target must be an ordinary file: $Path"
        }
    }
}
