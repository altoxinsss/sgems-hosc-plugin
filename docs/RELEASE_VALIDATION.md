# Release validation — 20 October 2026

The distribution provides prebuilt HOSC plugins for Windows x86 SGeMS 2.0,
Qt 3.3.5 and the VC7.1 runtime. Normal installation and example execution
require neither a compiler nor a separate Python installation.

## Package identity

| Artifact | SHA-256 |
|---|---|
| HOSC2008.dll | `5b98eba21ee4e46307f712a3160ad9bb92b5b6fe27b4fe95afe0e3868eda572d` |
| HOSCWindow.dll | `dbf265caabec37f87938c0b9dc368c1b3d80f2f7e5f345ec732c9ced087d6ba9` |

`bin/win32/manifest.json` also records both UI checksums and the tested host
fingerprints. `build_source_manifest.json` records the compiled source inputs.
The source was built in a separate directory containing spaces with the
VC7.1 compiler (13.10.3077); the numerical and Qt3 widget tests passed.

## Completed checks

| Boundary | Result |
|---|---|
| Installation | Preflight, first installation and repeated installation passed in a separate host directory. |
| Recovery | Fresh-install removal, upgrade restoration, corrupted-package rejection, modified-file protection and locked-file transaction rollback passed in isolated installation tests. |
| Plugin binary contract | Both DLLs are x86, import the VC7.1/Qt3 runtime and export `plugin_init`. |
| Core and input feedback | Numerical core checks and hidden Qt3 checks for field diagnostics, selected tab, text/combo highlighting and reset passed. |
| Host execution | Template and window calculations, exports and invalid-input recovery passed in one isolated SGeMS process, with exit code zero. |
| Supplied parameter files | All four `examples/*.par` files executed without changing their parameter values; 54 exported output nodes were checked. |
| Standalone source | Release build with MSVC 19.44 and CMake 3.31.6 passed; CTest passed 5/5 tests. |
| Published scientific evidence | 234 copied-file checksums, two derived CSVs and 11 figure/table contracts passed verification. |
| Download and Git checkout | ZIP extraction, temporary Git add/checkout with `core.autocrlf=true`, and execution from a path containing Unicode and spaces passed; package and research-evidence hashes remained valid. |

## Numerical examples

| Supplied parameter file | Calculation | Nodes | Maximum absolute raw-value error |
|---|---|---:|---:|
| `template_3d.par` | Fourth-order joint cumulant, AVL | 27 | 5.83e-12 |
| `window_uniform.par` | Fourth-order uniform cumulant, Bell | 9 | 2.43e-12 |
| `window_coordinate_gaussian.par` | Sixth-order spatially weighted cumulant | 9 | 2.10e-14 |
| `window_soft_high_porosity.par` | Sixth-order value-weighted cumulant | 9 | 1.13e-13 |

The template reference was calculated independently from centered fourth joint
moments and the three covariance products. Window results were checked against
`expected_window.csv`, including variance, standardized cumulant, ESS, support
and validity flags. Exported SGeMS properties use binary32 precision; audit
CSV values retain double precision. The reference tolerances are recorded in
each generated `RESULT.json`.

`RunExamples.cmd` reproduces these checks on the reviewer's host. It writes
the full parameter hashes, output checks and logs to a new result directory.
The two deliberately invalid calls produce four field diagnostics before
the valid cases run successfully. Interactive appearance of the complete
SGeMS application is distinct from the hidden Qt3 widget checks.

## External host and scientific-run scope

The official `sgems_2.0_installer.exe` was downloaded and its SHA-256 matched
the SourceForge metadata:
`9d1766826f0a24584e0877968789a6bce93636f0153afedbefc69ce908a7ab6f`.
The host used for execution was a separate copy of the working SGeMS 2.0
installation identified in `manifest.json`. The newly downloaded installer
was not installed during this verification. The installer checks binary
compatibility; the examples verify behavior on the selected installation.

The scientific result package was checked for integrity and provenance.
Full reservoir experiments were not repeated as part of deployment testing.
Its existing run reports describe the archived scientific results.

## Contents updated for this distribution

- Added the latest DLL/UI pairs, installation preflight, transactional backup
  and restoration, four loadable parameter files and executable examples.
- Replaced the default VC9 deployment route with the matching VC7.1 build.
  Removed the superseded VC9 build/deployment scripts and obsolete local
  `qstring.h` forwarding header; their previous copies were retained in an
  external development backup.
- Replaced development-machine log/temp paths with Windows temporary paths;
  removed an unused unsafe debug branch from template-result output helpers.
- Preserved the curated paper evidence; configured Git to retain its exact
  bytes and the checksummed UI/source files across different line-ending settings.
- Added portable source build, release staging and ZIP packaging tools.

The packaging command neither changes Git history nor publishes files.
