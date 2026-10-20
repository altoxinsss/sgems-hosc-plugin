param(
    [Parameter(Mandatory = $true)][string]$Vc71Root,
    [Parameter(Mandatory = $true)][string]$Qt3Dll,
    [Parameter(Mandatory = $true)][string]$VcVarsPath,
    [Parameter(Mandatory = $true)][string]$OutputDirectory
)

$ErrorActionPreference = "Stop"

function Invoke-QtNativeLogged([string]$Command, [string]$LogPath, [string]$Description) {
    # On Windows PowerShell 5.1, harmless native stderr becomes an error
    # record. Keep the log and use the process exit code as the result.
    $savedErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & $env:ComSpec /d /c $Command *> $LogPath
        $nativeExitCode = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
    if($nativeExitCode -ne 0) {
        throw "$Description failed with exit code $nativeExitCode. See $LogPath"
    }
}

$linker = Join-Path $Vc71Root "bin\link.exe"
$vcvars = $VcVarsPath
foreach($required in @($linker, $vcvars, $Qt3Dll)) {
    if(-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Qt3 import-library input was not found: $required"
    }
}
New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
$OutputDirectory = (Get-Item -LiteralPath $OutputDirectory).FullName
$Qt3Dll = (Get-Item -LiteralPath $Qt3Dll).FullName
$exportsLog = Join-Path $OutputDirectory "qt3_exports.log"
$linkLog = Join-Path $OutputDirectory "qt3_import_library.log"
$definitionPath = Join-Path $OutputDirectory "qt-mt335.def"
$libraryPath = Join-Path $OutputDirectory "qt-mt335.lib"
$manifestPath = Join-Path $OutputDirectory "qt3_import_manifest.json"

# Generate the import library from the selected host DLL so it matches that
# runtime, preserving both ordinals and DATA exports. In particular,
# qApp and QString::null must not be represented as callable import thunks.
$dumpCommand = 'call "{0}" >nul && "{1}" /dump /exports "{2}"' -f $vcvars, $linker, $Qt3Dll
Invoke-QtNativeLogged $dumpCommand $exportsLog "Qt3 export inspection"

$bytes = [System.IO.File]::ReadAllBytes($Qt3Dll)
$peOffset = [BitConverter]::ToInt32($bytes, 0x3c)
if([BitConverter]::ToUInt32($bytes, $peOffset) -ne 0x00004550 -or
   [BitConverter]::ToUInt16($bytes, $peOffset + 4) -ne 0x014c) {
    throw "The Qt3 runtime is not an x86 PE DLL: $Qt3Dll"
}
$sectionCount = [BitConverter]::ToUInt16($bytes, $peOffset + 6)
$optionalSize = [BitConverter]::ToUInt16($bytes, $peOffset + 20)
$sectionOffset = $peOffset + 24 + $optionalSize
$sections = @()
for($index = 0; $index -lt $sectionCount; ++$index) {
    $offset = $sectionOffset + 40 * $index
    $virtualSize = [BitConverter]::ToUInt32($bytes, $offset + 8)
    $rawSize = [BitConverter]::ToUInt32($bytes, $offset + 16)
    $sections += [PSCustomObject]@{
        rva = [BitConverter]::ToUInt32($bytes, $offset + 12)
        size = [Math]::Max($virtualSize, $rawSize)
        executable = (([BitConverter]::ToUInt32($bytes, $offset + 36) -band 0x20000000) -ne 0)
    }
}

$lines = New-Object 'System.Collections.Generic.List[string]'
$lines.Add(('LIBRARY "{0}"' -f [System.IO.Path]::GetFileName($Qt3Dll)))
$lines.Add("EXPORTS")
$exportCount = 0
$dataCount = 0
$qAppFound = $false
$declaredNames = 0
foreach($line in (Get-Content -LiteralPath $exportsLog)) {
    if($line -match '^\s*(\d+)\s+number of names\s*$') {
        $declaredNames = [int]$Matches[1]
    }
    if($line -notmatch '^\s*(\d+)\s+[0-9A-Fa-f]+\s+([0-9A-Fa-f]{8})\s+(\S+)\s*$') { continue }
    $ordinal = [int]$Matches[1]
    $rva = [Convert]::ToUInt32($Matches[2], 16)
    $symbol = $Matches[3]
    $section = @($sections | Where-Object { $rva -ge $_.rva -and $rva -lt ($_.rva + $_.size) })
    if($section.Count -ne 1) { throw "Could not classify Qt3 export $symbol at RVA $rva" }
    $isData = -not $section[0].executable
    $dataFlag = ""
    if($isData) { $dataFlag = " DATA"; ++$dataCount }
    if($symbol -eq '?qApp@@3PAVQApplication@@A') {
        if(-not $isData) { throw "The Qt3 qApp export was incorrectly classified as code." }
        $qAppFound = $true
    }
    $lines.Add(('    {0} @{1}{2}' -f $symbol, $ordinal, $dataFlag))
    ++$exportCount
}
if($exportCount -eq 0 -or $exportCount -ne $declaredNames -or -not $qAppFound) {
    throw "Incomplete Qt3 export table: parsed $exportCount of $declaredNames names; qApp DATA found=$qAppFound."
}
$lines | Set-Content -LiteralPath $definitionPath -Encoding ASCII
$libCommand = 'call "{0}" >nul && "{1}" /lib /nologo /machine:x86 /def:"{2}" /out:"{3}"' -f $vcvars, $linker, $definitionPath, $libraryPath
Invoke-QtNativeLogged $libCommand $linkLog "Qt3 import-library generation"
$manifest = [PSCustomObject]@{
    source_dll = $Qt3Dll
    source_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $Qt3Dll).Hash.ToLowerInvariant()
    machine = "x86"
    export_count = $exportCount
    data_export_count = $dataCount
    qApp_is_data = $qAppFound
    definition = $definitionPath
    library = $libraryPath
    library_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $libraryPath).Hash.ToLowerInvariant()
}
$manifest | ConvertTo-Json | Set-Content -LiteralPath $manifestPath -Encoding UTF8
Write-Output ([PSCustomObject]@{
    Library = $libraryPath
    Definition = $definitionPath
    Manifest = $manifestPath
    SourceDll = $Qt3Dll
    ExportCount = $exportCount
    DataExportCount = $dataCount
})
