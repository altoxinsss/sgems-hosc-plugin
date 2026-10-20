# Optional source build: VC7.1 x86

The prebuilt plugin package is the normal installation route for reviewers;
see the [repository README](../README.md). Compiling the SGeMS plugins is
optional. It requires separately obtained compiler, SDK, and Qt3 files that
are not part of this repository.

The supported binary boundary is the 32-bit SGeMS 2.0 host using Visual C++
7.1 (`MSVCR71.dll` and `MSVCP71.dll`) and `qt-mt335.dll`. The plugin exports
and the host exchange C++ objects. A newer compiler cannot be substituted
just because it can parse the source. In particular, **do not build with
Visual Studio 2008 / VC9 and install that DLL into the stock VC7.1 host**.
Old `.vcproj` files are historical material, not this build procedure.
The independent `standalone` implementation uses a modern C++17
toolchain and has its own [build instructions](../standalone/README.md).

## Required external inputs

| Input | Required contents |
|---|---|
| `Vc71Root` | VC7.1 compiler `bin/cl.exe` (version 13.10.x), linker `bin/link.exe` (7.10.x), C/C++ headers, x86 runtime import libraries, and a working x86 environment batch file |
| Windows SDK / Platform SDK | `Windows.h`, x86 `Kernel32.lib`, and the Windows headers/libraries selected by the VC7.1 environment batch file |
| `SGeMSSdkRoot` | Matching SGeMS 2.0 development tree: `GsTLAppli/`, `GsTL/`, `misc/`, the transitive dependency headers distributed with that SDK (including Boost), and `lib/` containing the matching x86 import libraries |
| `Qt3Root` | Qt **3.3.5** headers matching the host: either `QtKernel/`, `QtTools/`, `QtWidgets/`, or the usual `include/` directory |
| `Qt3Dll` | Absolute path to the target host's `qt-mt335.dll`; its directory must also contain that host's runtime dependencies, including `GsTLAppli_utils.dll` and the VC7.1 runtime DLLs |
| Shell | Windows PowerShell 5.1 or later; run in an ordinary user shell with write access to the build output directory |

The SGeMS import libraries required under `SGeMSSdkRoot/lib` are:

```text
coin2.lib                  GsTLAppli_actions.lib   GsTLAppli_appli.lib
GsTLAppli_extragui.lib      GsTLAppli_filters.lib   GsTLAppli_geostat.lib
GsTLAppli_grid.lib          GsTLAppli_gui.lib       GsTLAppli_math.lib
GsTLAppli_utils.lib         GsTLwidgets.lib         qwt.lib
simvoleon2.lib              soqt1.lib
```

Use the original development distribution corresponding to the installed
host, not headers from another SGeMS fork or import libraries from a newer
build. The compiler and its licensed runtime/import-library files must be
obtained separately through a permitted source. A SGeMS executable download
alone is not a development SDK. This repository does not download or
redistribute a compiler, SDK, Qt source tree, or third-party runtime DLLs.

The environment batch file must set `PATH`, `INCLUDE`, and `LIB` for **x86
VC7.1** and a compatible Windows SDK. The script detects these paths in order:
`Vc71Root/vcvars71.bat`, `Vc71Root/bin/vcvars32.bat`, and
`Vc71Root/vcvars32.bat`. Use `-VcVarsPath` if the environment file is elsewhere.
The build calls the compiler and linker under `Vc71Root/bin` explicitly.

## Build from any checkout directory

Replace the example external paths below with the locations on the build
machine. They are inputs, not prescribed installation directories. Paths
containing spaces are supported.

```powershell
Set-Location 'C:\work\sgems-hosc-plugin'

powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\build_v71_plugins.ps1 `
    -Vc71Root 'C:\toolchains\VC71' `
    -SGeMSSdkRoot 'C:\SDKs\SGeMS-2.0' `
    -Qt3Root 'C:\SDKs\Qt-3.3.5' `
    -Qt3Dll 'C:\Apps\SGeMS-2.0\qt-mt335.dll' `
    -StageDir '.\staging\reviewer_build'
```

For a nonstandard environment file, add:

```powershell
    -VcVarsPath 'C:\toolchains\VC71\configure-x86.bat'
```

`-LegacyRoot` remains an alias for `-SGeMSSdkRoot` for earlier local commands.
No source path or output directory needs the developer's original drive or
directory layout. Keep the checkout unchanged during the build: the script
compares source manifests before and after compiling.

## What a successful build verifies

The build stops on a failed command or a missing prerequisite. It:

1. Checks the selected compiler and linker versions and required SDK files.
2. Runs the template and window static interface-contract checks.
3. Generates a matching Qt3 import library from the selected host DLL,
   retaining export ordinals and marking data symbols such as `qApp` as
   `DATA`. The resulting library is generated in the staging directory.
4. Compiles and executes the shared weighted-statistics numerical checks.
5. Compiles and executes the Qt3 input-feedback checks in a separate hidden
   test process (field diagnostics, tab selection, highlighting, and reset).
6. Builds both plugin DLLs using x86 Release `/MD`, `/GR`, and `/EHsc`.
   `/GR` is required for SGeMS factory-object casts across the DLL boundary.
7. Parses the two Qt3 UI XML files, inspects the DLL exports/imports, verifies
   `plugin_init`, x86/VC7.1 linkage and the Qt3 runtime import, and confirms
   that the checked source manifest did not change during compilation.

The small native test source retains its historical name
`tests/vc9_core_checks.cpp`; in this workflow it is compiled and executed by
VC7.1. Static contract checks inspect source/UI structure; their success
alone does not establish numerical or host-runtime behavior.

## Outputs and installation

The chosen staging directory contains:

| Output | Purpose |
|---|---|
| `HOSC2008.dll` + `HOSC_TEMPLATE.ui` | Template plugin and matching interface |
| `HOSCWindow.dll` + `HOSC_WINDOW.ui` | Moving-window plugin and matching interface |
| `*.pdb` | Debug symbols for diagnosis; unnecessary for normal use |
| `vc71_core_checks.exe`, `vc71_validation_feedback_checks.exe` | Test programs run during the build |
| `logs/` | Per-source compile, link, test, PE-header, export, and import evidence |
| `source_manifest_prebuild.json`, `source_manifest_postbuild.json` | Relative source paths and SHA-256 hashes |
| `build_manifest.json`, `toolchain_report.md` | Artifact hashes, toolchain versions, and build summary |
| `qt3_import/` | Generated import library, export definition, and host-DLL hash |

Building does not install DLLs or change a running SGeMS instance. To use a
custom build, stage its matched DLL/UI files and regenerate the checksums:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\stage_release.ps1 `
    -StageDir '.\staging\reviewer_build' -SGeMSRoot 'C:\Apps\SGeMS-2.0'
.\Install.cmd 'C:\Apps\SGeMS-2.0'
.\RunExamples.cmd 'C:\Apps\SGeMS-2.0'
```

Close the selected SGeMS host before installation. Staging rejects a source
or artifact mismatch and records the compiled source hashes in
`bin/win32/build_source_manifest.json`. Installation creates a recoverable
backup. Use the DLL and UI from the **same build**. Additional fixture checks
are in [MANUAL_TESTING.md](MANUAL_TESTING.md).

Maintainers can create a ZIP of the current working files with
`powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\package_release.ps1`.
This command requires Git and includes new files while excluding removed
legacy trees and ignored build products. It does not commit, push, or publish.

The build-time checks do not run a full SGeMS workflow. Host load, calculation,
repeated execution, and save/reload checks remain separate. Source manifests
bind a build to its inputs; timestamped DLL/PDB binaries are not expected to
be byte-for-byte identical across independent builds.

## Troubleshooting

| Failure | Action |
|---|---|
| Compiler/environment script is missing | Set `Vc71Root` and, if needed, `VcVarsPath` to the actual VC7.1 installation. |
| `Windows.h` or `Kernel32.lib` cannot be found | Correct the Windows SDK paths in the external environment batch file. |
| A `GsTLAppli`, `GsTL`, or Boost header is missing | Supply the complete matching development SDK and its dependency headers. |
| Qt3 export parsing fails | Verify that `Qt3Dll` is the x86 `qt-mt335.dll` from the supported host, rather than Qt4/Qt5/Qt6 or a renamed DLL. |
| Qt3 hidden test cannot load a DLL | Check the host runtime directory and VC7.1 runtime availability; read `logs/vc71_validation_feedback_run.log`. |
| Link or runtime-library check fails | Read the named log; confirm all import libraries and `/MD` runtime libraries are x86 VC7.1 compatible. |
| Source manifest changed during build | Stop other source-editing tasks and rebuild into a fresh staging directory. |

Keep the error log and `toolchain_report.md` when reporting a build issue.
