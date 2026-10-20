param(
    [Parameter(Mandatory = $true)][string]$Vc71Root,
    [Parameter(Mandatory = $true)][Alias("LegacyRoot")][string]$SGeMSSdkRoot,
    [Parameter(Mandatory = $true)][string]$Qt3Root,
    [Parameter(Mandatory = $true)][string]$Qt3Dll,
    [string]$VcVarsPath = "",
    [string]$StageDir = ""
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $PSScriptRoot

function Require-Path([string]$Path, [string]$Description) {
    if(-not (Test-Path -LiteralPath $Path)) {
        throw "$Description was not found: $Path"
    }
}

function Invoke-CmdLogged([string]$Command, [string]$LogPath, [string]$Description) {
    # Windows PowerShell 5.1 wraps native stderr as NativeCommandError even
    # when the executable succeeds (the VC7.1 banner is written to stderr).
    # Preserve both streams, then decide success from the native exit code.
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

function Get-Sha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Add-ManifestEntry([System.Collections.Generic.List[object]]$Entries, [string]$Path) {
    $item = Get-Item -LiteralPath $Path
    $Entries.Add([PSCustomObject]@{
        path = $item.FullName
        sha256 = Get-Sha256 $item.FullName
        bytes = $item.Length
    })
}

function Get-SourceManifest([string]$RepositoryRoot) {
    $roots = @("compat", "shared", "MoveWindow\HOSCWindow", "Template\HOSC2008", "ui", "tests", "tools", "docs") |
        ForEach-Object { Join-Path $RepositoryRoot $_ }
    $extensions = @(".c", ".cc", ".cpp", ".h", ".vcproj", ".sln", ".ui", ".ps1", ".md", ".sgems")
    $entries = @()
    foreach($root in $roots) {
        if(-not (Test-Path -LiteralPath $root)) { continue }
        Get-ChildItem -LiteralPath $root -File -Recurse |
            Where-Object { $extensions -contains $_.Extension.ToLowerInvariant() } |
            ForEach-Object {
                $entries += [PSCustomObject]@{
                    relative_path = $_.FullName.Substring($RepositoryRoot.Length).TrimStart('\')
                    bytes = $_.Length
                    sha256 = Get-Sha256 $_.FullName
                }
            }
    }
    return @($entries | Sort-Object relative_path)
}

function Get-CommandArgument([string]$Argument) {
    return '"{0}"' -f $Argument
}

function Build-Plugin(
    [string]$Name,
    [string]$SourceDirectory,
    [string]$OutputDirectory,
    [string]$LogsDirectory,
    [string]$VcVars,
    [string]$Compiler,
    [string]$Linker,
    [string]$IncludeFlags,
    [string]$LibraryDirectory,
    [string]$Dependencies
) {
    New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
    # S-GeMS creates plugins through Named_interface and dynamically casts them
    # back to Geostat_algo across the DLL boundary.  VC7.1 otherwise defaults to
    # /GR-, which leaves the factory object constructible but prevents the host
    # from reaching initialize().
    $compileFlags = "/nologo /EHsc /GR /MD /O2 /W3 /Zi /DWIN32 /DNDEBUG /D_CONSOLE /DQT_DLL /DQT_SHARED /DQT_CORE_LIB /DQT_GUI_LIB /D_WINDOWS $IncludeFlags"
    $objects = @()
    $sources = @(
        "AVL.cpp",
        "CumulantCal.cpp",
        "FunctionCall.cpp",
        "Initialization.cpp",
        "PCA.cpp",
        "Plugin.cpp",
        "SearchNode.cpp",
        "Tools.cpp"
    )
    foreach($sourceName in $sources) {
        $source = Join-Path $SourceDirectory $sourceName
        Require-Path $source "$Name source file $sourceName"
        $object = Join-Path $OutputDirectory ([System.IO.Path]::GetFileNameWithoutExtension($sourceName) + ".obj")
        $log = Join-Path $LogsDirectory ("{0}_{1}.log" -f $Name, [System.IO.Path]::GetFileNameWithoutExtension($sourceName))
        $command = 'call {0} >nul && cd /d {1} && {2} {3} /Fo{4} /c {5}' -f `
            (Get-CommandArgument $VcVars),
            (Get-CommandArgument $OutputDirectory),
            (Get-CommandArgument $Compiler),
            $compileFlags,
            (Get-CommandArgument $object),
            (Get-CommandArgument $source)
        Invoke-CmdLogged $command $log "$Name compile $sourceName"
        $objects += $object
    }

    $sharedSource = Join-Path $repoRoot "shared\WeightedMomentCore.cpp"
    Require-Path $sharedSource "Shared weighted-moment source"
    $sharedObject = Join-Path $OutputDirectory "WeightedMomentCore.obj"
    $sharedLog = Join-Path $LogsDirectory ("{0}_WeightedMomentCore.log" -f $Name)
    $sharedCommand = 'call {0} >nul && cd /d {1} && {2} {3} /Fo{4} /c {5}' -f `
        (Get-CommandArgument $VcVars),
        (Get-CommandArgument $OutputDirectory),
        (Get-CommandArgument $Compiler),
        $compileFlags,
        (Get-CommandArgument $sharedObject),
        (Get-CommandArgument $sharedSource)
    Invoke-CmdLogged $sharedCommand $sharedLog "$Name compile shared weighted-moment core"
    $objects += $sharedObject

    $dll = Join-Path $OutputDirectory ($Name + ".dll")
    $pdb = Join-Path $OutputDirectory ($Name + ".pdb")
    $linkLog = Join-Path $LogsDirectory ($Name + "_link.log")
    $objectArgs = ($objects | ForEach-Object { Get-CommandArgument $_ }) -join " "
    $linkCommand = 'call {0} >nul && cd /d {1} && {2} /nologo /DLL /MACHINE:X86 /DEBUG /OUT:{3} /PDB:{4} {5} /LIBPATH:{6} {7}' -f `
        (Get-CommandArgument $VcVars),
        (Get-CommandArgument $OutputDirectory),
        (Get-CommandArgument $Linker),
        (Get-CommandArgument $dll),
        (Get-CommandArgument $pdb),
        $objectArgs,
        (Get-CommandArgument $LibraryDirectory),
        $Dependencies
    Invoke-CmdLogged $linkCommand $linkLog "$Name DLL link"
    Require-Path $dll "$Name output DLL"
    Require-Path $pdb "$Name output PDB"
    return [PSCustomObject]@{ Dll = $dll; Pdb = $pdb }
}

foreach($inputDirectory in @($Vc71Root, $SGeMSSdkRoot, $Qt3Root)) {
    Require-Path $inputDirectory "Build input directory"
}
$Vc71Root = (Get-Item -LiteralPath $Vc71Root).FullName
$SGeMSSdkRoot = (Get-Item -LiteralPath $SGeMSSdkRoot).FullName
$Qt3Root = (Get-Item -LiteralPath $Qt3Root).FullName
Require-Path $Qt3Dll "Host Qt3 runtime DLL"
$Qt3Dll = (Get-Item -LiteralPath $Qt3Dll).FullName
if([string]::IsNullOrWhiteSpace($VcVarsPath)) {
    foreach($candidate in @("vcvars71.bat", "bin\vcvars32.bat", "vcvars32.bat")) {
        $candidatePath = Join-Path $Vc71Root $candidate
        if(Test-Path -LiteralPath $candidatePath -PathType Leaf) {
            $VcVarsPath = $candidatePath
            break
        }
    }
    if([string]::IsNullOrWhiteSpace($VcVarsPath)) {
        throw "No VC7.1 environment script was found. Supply -VcVarsPath pointing to an x86 VC7.1 environment .bat file."
    }
}
Require-Path $VcVarsPath "VC7.1 environment script"
$vcvars = (Get-Item -LiteralPath $VcVarsPath).FullName
$compiler = Join-Path $Vc71Root "bin\cl.exe"
$linker = Join-Path $Vc71Root "bin\link.exe"
$legacyLibDirectory = Join-Path $SGeMSSdkRoot "lib"

Require-Path $vcvars "VC7.1 environment script"
Require-Path $compiler "VC7.1 compiler"
Require-Path $linker "VC7.1 linker"
Require-Path $legacyLibDirectory "Matching SGeMS SDK library directory"
Require-Path (Join-Path $SGeMSSdkRoot "GsTLAppli\geostat\geostat_algo.h") "SGeMS SDK header"
Require-Path (Join-Path $SGeMSSdkRoot "GsTL\math\math_functions.h") "GsTL dependency header"
Require-Path (Join-Path $SGeMSSdkRoot "misc") "SGeMS misc include directory"
$qt3IncludeDirectories = @()
if(Test-Path -LiteralPath (Join-Path $Qt3Root "QtKernel\qwidget.h")) {
    $qt3IncludeDirectories = @("QtKernel", "QtTools", "QtWidgets") | ForEach-Object { Join-Path $Qt3Root $_ }
    $qtGlobalHeader = Join-Path $Qt3Root "QtTools\qglobal.h"
    Require-Path $qtGlobalHeader "Qt3 global header"
    Require-Path (Join-Path $Qt3Root "QtWidgets\qtextedit.h") "Qt3 text-edit header"
}
elseif(Test-Path -LiteralPath (Join-Path $Qt3Root "include\qwidget.h")) {
    $qt3IncludeDirectories = @(Join-Path $Qt3Root "include")
    $qtGlobalHeader = Join-Path $Qt3Root "include\qglobal.h"
    Require-Path $qtGlobalHeader "Qt3 global header"
    Require-Path (Join-Path $Qt3Root "include\qtextedit.h") "Qt3 text-edit header"
}
else {
    throw "Qt3Root must contain QtKernel/QtTools/QtWidgets or an include directory with the Qt 3.3.5 headers."
}
if((Get-Content -LiteralPath $qtGlobalHeader -Raw) -notmatch '(?m)^\s*#\s*define\s+QT_VERSION\s+0x030305\b') {
    throw "The selected Qt headers do not declare Qt 3.3.5: $qtGlobalHeader"
}
if([System.IO.Path]::GetFileName($Qt3Dll) -ine "qt-mt335.dll") {
    throw "Qt3Dll must identify the target host's qt-mt335.dll."
}
Require-Path $Qt3Dll "Host Qt3 runtime DLL"
Require-Path (Join-Path $repoRoot "tools\prepare_qt3_import_library.ps1") "Qt3 import-library generator"
foreach($library in @("coin2.lib", "GsTLAppli_actions.lib", "GsTLAppli_appli.lib", "GsTLAppli_extragui.lib", "GsTLAppli_filters.lib", "GsTLAppli_geostat.lib", "GsTLAppli_grid.lib", "GsTLAppli_gui.lib", "GsTLAppli_math.lib", "GsTLAppli_utils.lib", "GsTLwidgets.lib", "qwt.lib", "simvoleon2.lib", "soqt1.lib")) {
    Require-Path (Join-Path $legacyLibDirectory $library) "Legacy dependency $library"
}

$compilerVersion = (Get-Item -LiteralPath $compiler).VersionInfo.FileVersion
if($compilerVersion -notmatch "^13\.10\.") {
    throw "The selected compiler is not VC7.1: $compiler"
}
$linkerVersion = (Get-Item -LiteralPath $linker).VersionInfo.FileVersion
if($linkerVersion -notmatch "^7\.10\.") {
    throw "The selected linker is not VC7.1: $linker"
}

if([string]::IsNullOrEmpty($StageDir)) {
    $StageDir = Join-Path $repoRoot ("staging\vc71_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
}
New-Item -ItemType Directory -Force -Path $StageDir | Out-Null
$StageDir = (Get-Item -LiteralPath $StageDir).FullName
$logsDirectory = Join-Path $StageDir "logs"
New-Item -ItemType Directory -Force -Path $logsDirectory | Out-Null

$compilerProbeLog = Join-Path $logsDirectory "vc71_compiler_probe.log"
# cl /Bv intentionally returns D2003 when no translation unit is supplied.
# Preserve its version output while normalizing that expected diagnostic to a
# successful probe command.
Invoke-CmdLogged ('call {0} >nul && {1} /Bv & exit /b 0' -f (Get-CommandArgument $vcvars), (Get-CommandArgument $compiler)) $compilerProbeLog "VC7.1 compiler probe"
if((Get-Content -LiteralPath $compilerProbeLog -Raw) -notmatch "13\.10\.") {
    throw "The compiler probe did not report a VC7.1 (13.10.x) compiler."
}

& (Join-Path $repoRoot "tools\verify_template_contract.ps1") -RepoRoot $repoRoot
& (Join-Path $repoRoot "tools\verify_moving_window_contract.ps1") -RepoRoot $repoRoot

$preBuildManifest = Get-SourceManifest $repoRoot
$preBuildManifestPath = Join-Path $StageDir "source_manifest_prebuild.json"
$preBuildManifest | ConvertTo-Json | Set-Content -LiteralPath $preBuildManifestPath -Encoding UTF8

$qt3Import = & (Join-Path $repoRoot "tools\prepare_qt3_import_library.ps1") `
    -Vc71Root $Vc71Root -VcVarsPath $vcvars -Qt3Dll $Qt3Dll -OutputDirectory (Join-Path $StageDir "qt3_import")

$includeDirectories = @(
    (Join-Path $repoRoot "compat\vc71"),
    $SGeMSSdkRoot,
    (Join-Path $SGeMSSdkRoot "misc"),
    (Join-Path $repoRoot "Template"),
    (Join-Path $repoRoot "Template\HOSC2008"),
    (Join-Path $repoRoot "MoveWindow\HOSCWindow"),
    (Join-Path $repoRoot "shared")
) + $qt3IncludeDirectories
$includeFlags = (($includeDirectories | ForEach-Object { '/I "{0}"' -f $_ }) -join " ")
$dependencies = "coin2.lib GsTLAppli_actions.lib GsTLAppli_appli.lib GsTLAppli_extragui.lib GsTLAppli_filters.lib GsTLAppli_geostat.lib GsTLAppli_grid.lib GsTLAppli_gui.lib GsTLAppli_math.lib GsTLAppli_utils.lib GsTLwidgets.lib qwt.lib simvoleon2.lib soqt1.lib"
$dependencies += " " + (Get-CommandArgument $qt3Import.Library)

$coreDirectory = Join-Path $StageDir "core_checks"
New-Item -ItemType Directory -Force -Path $coreDirectory | Out-Null
$coreCompileLog = Join-Path $logsDirectory "vc71_core_checks_compile.log"
$coreLinkLog = Join-Path $logsDirectory "vc71_core_checks_link.log"
$coreRunLog = Join-Path $logsDirectory "vc71_core_checks_run.log"
$coreSource = Join-Path $repoRoot "tests\vc9_core_checks.cpp"
$sharedSource = Join-Path $repoRoot "shared\WeightedMomentCore.cpp"
$coreCompileCommand = 'call {0} >nul && cd /d {1} && {2} /nologo /EHsc /MD /O2 /W4 /I {3} /c {4} {5}' -f `
    (Get-CommandArgument $vcvars),
    (Get-CommandArgument $coreDirectory),
    (Get-CommandArgument $compiler),
    (Get-CommandArgument (Join-Path $repoRoot "shared")),
    (Get-CommandArgument $sharedSource),
    (Get-CommandArgument $coreSource)
Invoke-CmdLogged $coreCompileCommand $coreCompileLog "VC7.1 weighted-core compile"
$coreExecutable = Join-Path $StageDir "vc71_core_checks.exe"
$coreLinkCommand = 'call {0} >nul && cd /d {1} && {2} /nologo WeightedMomentCore.obj vc9_core_checks.obj /OUT:{3}' -f `
    (Get-CommandArgument $vcvars),
    (Get-CommandArgument $coreDirectory),
    (Get-CommandArgument $linker),
    (Get-CommandArgument $coreExecutable)
Invoke-CmdLogged $coreLinkCommand $coreLinkLog "VC7.1 weighted-core link"
Invoke-CmdLogged ('call {0} >nul && {1}' -f (Get-CommandArgument $vcvars), (Get-CommandArgument $coreExecutable)) $coreRunLog "VC7.1 weighted-core checks"

# Exercise real Qt3 widgets in a hidden form. This separate process never
# touches an existing S-GeMS window and does not display modal message boxes.
$validationDirectory = Join-Path $StageDir "validation_checks"
New-Item -ItemType Directory -Force -Path $validationDirectory | Out-Null
$validationSource = Join-Path $repoRoot "tests\validation_feedback_checks.cpp"
$validationExecutable = Join-Path $StageDir "vc71_validation_feedback_checks.exe"
$validationCompileLog = Join-Path $logsDirectory "vc71_validation_feedback_compile.log"
$validationLinkLog = Join-Path $logsDirectory "vc71_validation_feedback_link.log"
$validationRunLog = Join-Path $logsDirectory "vc71_validation_feedback_run.log"
Require-Path $validationSource "Qt3 validation-feedback test source"
$validationCompileCommand = 'call {0} >nul && cd /d {1} && {2} /nologo /EHsc /GR /MD /O2 /W4 /DWIN32 /DQT_DLL {3} /c {4}' -f `
    (Get-CommandArgument $vcvars), (Get-CommandArgument $validationDirectory), `
    (Get-CommandArgument $compiler), $includeFlags, (Get-CommandArgument $validationSource)
Invoke-CmdLogged $validationCompileCommand $validationCompileLog "VC7.1 Qt3 validation-feedback compile"
$validationLinkCommand = 'call {0} >nul && cd /d {1} && {2} /nologo validation_feedback_checks.obj /OUT:{3} /LIBPATH:{4} GsTLAppli_utils.lib {5}' -f `
    (Get-CommandArgument $vcvars), (Get-CommandArgument $validationDirectory), `
    (Get-CommandArgument $linker), (Get-CommandArgument $validationExecutable), `
    (Get-CommandArgument $legacyLibDirectory), (Get-CommandArgument $qt3Import.Library)
Invoke-CmdLogged $validationLinkCommand $validationLinkLog "VC7.1 Qt3 validation-feedback link"
$validationSavedPath = $env:PATH
$validationSavedErrorActionPreference = $ErrorActionPreference
try {
    # vcvars71.bat retains a saved PATH if called by an already configured
    # shell; invoke the linked test directly so the runtime directory is first.
    $env:PATH = (Split-Path -Parent $Qt3Dll) + ";" + (Join-Path $Vc71Root "bin") + ";" + $validationSavedPath
    $ErrorActionPreference = "Continue"
    & $validationExecutable *> $validationRunLog
    $validationExitCode = $LASTEXITCODE
    $ErrorActionPreference = $validationSavedErrorActionPreference
    if($validationExitCode -ne 0) {
        throw "VC7.1 Qt3 validation-feedback checks failed with exit code $validationExitCode. See $validationRunLog"
    }
}
finally {
    $env:PATH = $validationSavedPath
    $ErrorActionPreference = $validationSavedErrorActionPreference
}

$template = Build-Plugin -Name "HOSC2008" `
    -SourceDirectory (Join-Path $repoRoot "Template\HOSC2008") `
    -OutputDirectory (Join-Path $StageDir "obj\HOSC2008") `
    -LogsDirectory $logsDirectory `
    -VcVars $vcvars -Compiler $compiler -Linker $linker `
    -IncludeFlags $includeFlags -LibraryDirectory $legacyLibDirectory -Dependencies $dependencies
$window = Build-Plugin -Name "HOSCWindow" `
    -SourceDirectory (Join-Path $repoRoot "MoveWindow\HOSCWindow") `
    -OutputDirectory (Join-Path $StageDir "obj\HOSCWindow") `
    -LogsDirectory $logsDirectory `
    -VcVars $vcvars -Compiler $compiler -Linker $linker `
    -IncludeFlags $includeFlags -LibraryDirectory $legacyLibDirectory -Dependencies $dependencies

Copy-Item -LiteralPath $template.Dll -Destination (Join-Path $StageDir "HOSC2008.dll") -Force
Copy-Item -LiteralPath $template.Pdb -Destination (Join-Path $StageDir "HOSC2008.pdb") -Force
Copy-Item -LiteralPath $window.Dll -Destination (Join-Path $StageDir "HOSCWindow.dll") -Force
Copy-Item -LiteralPath $window.Pdb -Destination (Join-Path $StageDir "HOSCWindow.pdb") -Force

foreach($uiName in @("HOSC_TEMPLATE.ui", "HOSC_WINDOW.ui")) {
    $source = Join-Path $repoRoot ("ui\" + $uiName)
    Require-Path $source "UI source $uiName"
    $ui = New-Object System.Xml.XmlDocument
    try {
        $ui.Load($source)
    }
    catch {
        throw "$uiName is not well-formed Qt3 UI XML: $($_.Exception.Message)"
    }
    if($ui.DocumentElement.Name -ne "UI" -or $ui.DocumentElement.GetAttribute("version") -ne "3.3") {
        throw "$uiName is not a Qt3.3 UI document."
    }
    Copy-Item -LiteralPath $source -Destination (Join-Path $StageDir $uiName) -Force
}

foreach($dllName in @("HOSC2008.dll", "HOSCWindow.dll")) {
    $dll = Join-Path $StageDir $dllName
    $headersLog = Join-Path $logsDirectory ($dllName + ".headers.log")
    $importsLog = Join-Path $logsDirectory ($dllName + ".imports.log")
    $exportsLog = Join-Path $logsDirectory ($dllName + ".exports.log")
    Invoke-CmdLogged ('call {0} >nul && {1} /dump /headers {2}' -f (Get-CommandArgument $vcvars), (Get-CommandArgument $linker), (Get-CommandArgument $dll)) $headersLog "$dllName header inspection"
    Invoke-CmdLogged ('call {0} >nul && {1} /dump /imports {2}' -f (Get-CommandArgument $vcvars), (Get-CommandArgument $linker), (Get-CommandArgument $dll)) $importsLog "$dllName import inspection"
    Invoke-CmdLogged ('call {0} >nul && {1} /dump /exports {2}' -f (Get-CommandArgument $vcvars), (Get-CommandArgument $linker), (Get-CommandArgument $dll)) $exportsLog "$dllName export inspection"
    $headers = Get-Content -LiteralPath $headersLog -Raw
    $imports = Get-Content -LiteralPath $importsLog -Raw
    $exports = Get-Content -LiteralPath $exportsLog -Raw
    if($headers -notmatch "machine \(x86\)" -or $headers -notmatch "7\.10") {
        throw "$dllName is not a VC7.1 x86 binary."
    }
    if($imports -notmatch "MSVCR71\.dll" -or $imports -notmatch "MSVCP71\.dll" -or $imports -match "MSVCR90|MSVCP90|MSVCR1[4-9]|MSVCP1[4-9]") {
        throw "$dllName does not have the required VC7.1 runtime dependencies."
    }
    if($exports -notmatch "plugin_init") {
        throw "$dllName does not export plugin_init."
    }
    if($imports -notmatch "qt-mt335\.dll") {
        throw "$dllName does not import the host Qt3 runtime required for validation feedback."
    }
}

$postBuildManifest = Get-SourceManifest $repoRoot
$postBuildManifestPath = Join-Path $StageDir "source_manifest_postbuild.json"
$postBuildManifest | ConvertTo-Json | Set-Content -LiteralPath $postBuildManifestPath -Encoding UTF8
if(($preBuildManifest | ConvertTo-Json -Compress) -ne ($postBuildManifest | ConvertTo-Json -Compress)) {
    throw "Source files changed during the VC7.1 build; staged binaries are not bound to one frozen source manifest."
}

$toolchainReport = Join-Path $StageDir "toolchain_report.md"
@(
    "# VC7.1 Plugin Build Report",
    "",
    "- Compiler: $compilerVersion",
    "- Linker: $((Get-Item -LiteralPath $linker).VersionInfo.FileVersion)",
    "- Runtime contract: x86, Release, /MD, MSVCR71.dll, MSVCP71.dll",
    "- SGeMS SDK headers and import libraries: $SGeMSSdkRoot",
    "- Host Qt3 runtime: $Qt3Dll",
    "- Qt3 import library: $($qt3Import.ExportCount) named exports, $($qt3Import.DataExportCount) DATA exports",
    "- Hidden Qt3 validation-feedback checks: passed (parse, tab selection, editable and combo highlighting, reset)",
    "- Staging directory: $StageDir",
    "- Source manifest SHA-256: $(Get-Sha256 $preBuildManifestPath)",
    "",
    "Both plugin DLLs were rebuilt with the same VC7.1 C++ runtime as the S-GeMS host boundary.",
    "Deployment is intentionally separate from this build script."
) | Set-Content -LiteralPath $toolchainReport -Encoding UTF8

$entries = New-Object 'System.Collections.Generic.List[object]'
foreach($path in @(
    $coreExecutable,
    $validationExecutable,
    $qt3Import.Library,
    $qt3Import.Definition,
    $qt3Import.Manifest,
    (Join-Path $StageDir "HOSC2008.dll"),
    (Join-Path $StageDir "HOSC2008.pdb"),
    (Join-Path $StageDir "HOSCWindow.dll"),
    (Join-Path $StageDir "HOSCWindow.pdb"),
    (Join-Path $StageDir "HOSC_TEMPLATE.ui"),
    (Join-Path $StageDir "HOSC_WINDOW.ui"),
    $preBuildManifestPath,
    $postBuildManifestPath,
    $toolchainReport
)) {
    Add-ManifestEntry $entries $path
}
$entries | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $StageDir "build_manifest.json") -Encoding UTF8

Write-Host "VC7.1 plugin build completed successfully. Staged artifacts: $StageDir"
