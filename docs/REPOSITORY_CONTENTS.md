# Repository contents

## Installable package and examples

| Path | Purpose |
|---|---|
| `Install.cmd` | Interactive or command-line installation entry point |
| `RunExamples.cmd` | Run and check small examples in an isolated SGeMS environment |
| `bin/win32/` | Prebuilt x86 VC7.1 plugin DLLs, matching interfaces, and package manifest |
| `tools/install_plugin.ps1` | Host/package preflight, backup, and deployment |
| `tools/uninstall_plugin.ps1` | Restore an installation backup |
| `tools/run_examples.ps1` | Example execution, exported outputs, and numerical validation |
| `tools/stage_release.ps1` | Bind rebuilt DLL/UI artifacts to the compiled source hashes |
| `tools/package_release.ps1` | Create a ZIP and checksum inventory of the current working files |
| `examples/*.par` | Loadable SGeMS parameter files for exact templates and uniform/spatial/value-weighted windows |
| `tests/fixtures/manual/selftest_3d_cartesian.sgems` | Native `5 × 5 × 3` Cartesian grid named `test`, property `phi` |
| `tests/fixtures/manual/*.gslib` | Additional synthetic properties for manual integration tests |
| `tests/fixtures/manual/expected_*.csv` | Numerical reference values for the documented fixture cases |

No compiler or Python installation is required to install the prebuilt package
and run its SGeMS examples. The matching SGeMS application is an external
prerequisite; it is not copied into this repository.

## Source code

| Path | Purpose |
|---|---|
| `Template/HOSC2008/` | Template plugin implementation; the directory name is retained from the original project |
| `MoveWindow/HOSCWindow/` | Moving-window plugin implementation |
| `shared/` | Normalized weighted-statistics core compiled into both plugins |
| `compat/vc71/` | Compiler compatibility support for the VC7.1 build |
| `ui/` | Source Qt3 interfaces for `HOSC_TEMPLATE` and `HOSC_WINDOW` |
| `tools/build_v71_plugins.ps1` | Optional plugin build using an externally supplied matching development environment |
| `standalone/` | Independent C++17 command-line calculations, native checks, and Python experiment workflows |

The prebuilt plugins use VC7.1 to match the supported SGeMS host. Retained
Visual Studio project filenames do not change the binary compatibility
requirement. See [the build guide](BUILD_FROM_SOURCE.md) before rebuilding.

## Scientific data and reproduction

| Path | Purpose |
|---|---|
| `Paper-data/current_results/2026-10-12/` | Compact results and provenance for Sections 4.3–4.6 and the conditional-shape supplement |
| `Paper-data/tools/` | Evidence integrity checks and figure-input export |
| `Paper-data/raw_data/README.txt` | Stanford VI-E attribution and source-data preparation instructions |
| `standalone/experiments/` | Orders 4–8 benchmarks, directional templates, direct window statistics, controlled-moment checks, and supporting workflows |
| `standalone/configs/weighting/` | Uniform, coordinate-Gaussian, and soft-high-porosity configuration examples |

The `Paper-data` records preserve the original scientific-run configurations,
reports, and hashes. Use the reproduction guide to resolve source-data and
output paths on a new machine. Large raw reservoir files and generated maps
are obtained or generated separately.

Some module directory names retain their original dates. In particular,
`standalone/experiments/r1_5_incremental_value_20260924/` supplies modules used
by the direct-statistics workflow, and `standalone/vendor/` contains active
numerical build dependencies.

## Documentation and verification

| Document | Purpose |
|---|---|
| [QUICKSTART.md](QUICKSTART.md) | Installation, worked interface settings, and expected outputs |
| [TROUBLESHOOTING.md](TROUBLESHOOTING.md) | Host, parameter, and output diagnosis |
| [NUMERICAL_CONTRACT.md](NUMERICAL_CONTRACT.md) | Statistical definitions, support, geometry, and audit rules |
| [PARAMETER_MAPPING.md](PARAMETER_MAPPING.md) | Interface and implementation parameter correspondence |
| [MANUAL_TESTING.md](MANUAL_TESTING.md) | Extended numerical, invalid-input, and persistence checks |
| [Tests overview](../tests/README.md) | Test layers and reference-fixture scope |
| [Release validation](RELEASE_VALIDATION.md) | Checks performed on this distribution and their scope |
| [Standalone reproduction](../standalone/experiments/README.md) | Full scientific-workflow requirements |

An example validation manifest records a particular execution on the selected
host. Core tests, static checks, and scientific evidence checks address separate
boundaries; their outputs identify what was tested. The small synthetic
examples are software checks rather than manuscript experimental results.

Dependencies and data retain their applicable notices; see
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
