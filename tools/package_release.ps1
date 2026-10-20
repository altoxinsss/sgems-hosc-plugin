#Requires -Version 5.1
[CmdletBinding()]
param([string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'host_common.ps1')
$package = Get-HoscPackage (Join-Path $repo 'bin\win32')
if (!$OutputDirectory) { $OutputDirectory = Join-Path $repo 'release-packages' }
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$OutputDirectory = (Get-Item -LiteralPath $OutputDirectory).FullName
$name = 'sgems-hosc-plugin-' + (Get-Date -Format 'yyyyMMdd_HHmmss')
$archive = Join-Path $OutputDirectory ($name + '.zip')
if (Test-Path -LiteralPath $archive) { throw "Archive already exists: $archive" }
Push-Location $repo
try {
    # Package the CURRENT working files, including new release artifacts and
    # excluding deleted legacy files. HEAD may predate the release preparation.
    $files = @(& git -c core.quotepath=false ls-files --cached --others --exclude-standard)
    if ($LASTEXITCODE -ne 0) { throw 'Git could not enumerate release files.' }
    $files = @($files | Sort-Object -Unique | Where-Object {
        $_ -notmatch '^(release-packages|reviewer-output|staging|evidence|output)/' -and
        (Test-Path -LiteralPath $_ -PathType Leaf)
    })
    foreach ($required in @('Install.cmd','RunExamples.cmd','README.md','bin/win32/manifest.json',
            'bin/win32/HOSC2008.dll','bin/win32/HOSCWindow.dll',
            'bin/win32/HOSC_TEMPLATE.ui','bin/win32/HOSC_WINDOW.ui')) {
        if ($files -notcontains $required) { throw "Required file is excluded from Git: $required" }
    }
    Add-Type -AssemblyName System.IO.Compression
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip = [IO.Compression.ZipFile]::Open($archive, [IO.Compression.ZipArchiveMode]::Create)
    $manifest = New-Object 'System.Collections.Generic.List[object]'
    try {
        foreach ($relative in $files) {
            $source = Join-Path $repo $relative
            $item = Get-Item -LiteralPath $source
            if ($item.Attributes -band [IO.FileAttributes]::ReparsePoint) {
                throw "Release files must be regular files: $relative"
            }
            if ($item.Length -ge 100MB) { throw "File exceeds GitHub's regular-file limit: $relative" }
            [IO.Compression.ZipFileExtensions]::CreateEntryFromFile($zip, $source,
                ('sgems-hosc-plugin/' + $relative.Replace('\','/')),
                [IO.Compression.CompressionLevel]::Optimal) | Out-Null
            $manifest.Add([PSCustomObject]@{
                path = $relative.Replace('\','/'); bytes = $item.Length
                sha256 = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash.ToLowerInvariant()
            })
        }
    }
    finally { $zip.Dispose() }
    $manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $OutputDirectory ($name + '.files.json')) -Encoding UTF8
    $hash = (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant()
    "$hash  $([IO.Path]::GetFileName($archive))" | Set-Content -LiteralPath ($archive + '.sha256') -Encoding ASCII
    Write-Host "Packaged $($files.Count) current files: $archive"
    Write-Host "SHA256: $hash"
}
finally { Pop-Location }
