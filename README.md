# HOSC plugins for SGeMS

Calculate high-order spatial-template and moving-window statistics on regular
Cartesian grids in SGeMS 2.0. This repository includes ready-to-install Windows
plugins, matching interfaces, small executable examples, source code, and the
paper's compact numerical evidence.

## Install and run

The prebuilt plugins target **32-bit SGeMS 2.0 with Qt 3.3.5 and the VC7.1
runtime** on Windows. A 64-bit Windows system can run this 32-bit application.
No compiler or Python installation is needed for this route.

1. Download this repository using **Code → Download ZIP** and extract it, or
   clone it with Git.
2. Install a compatible Windows SGeMS 2.0 distribution from the
   [SGeMS 2.0 Windows installer directory](https://sourceforge.net/projects/sgems/files/sgems/v2.0%20Windows%20installer/)
   (`sgems_2.0_installer.exe`), then close
   SGeMS before installing the plugins. The SGeMS application and its runtime
   dependencies are supplied by the SGeMS distribution.
3. Double-click `Install.cmd` and provide the folder containing `sgems.exe`.
   The installer checks the host and package, backs up replaced plugin files,
   and installs both DLL/UI pairs.
4. Run `RunExamples.cmd` with the same SGeMS folder to calculate and validate
   the supplied examples in a separate test directory.
5. Open SGeMS and select `HOSC_WINDOW` or `HOSC_TEMPLATE` from its algorithms.

For example, from a PowerShell window opened in this repository:

```powershell
.\Install.cmd "C:\SGeMS"
.\RunExamples.cmd "C:\SGeMS"
```

[Quick start](docs/QUICKSTART.md) provides the complete interface settings,
expected values, data-import instructions, and output descriptions.
[Interface guide](docs/UI_GUIDE.md) shows the current panels and explains
which controls apply to each calculation.
[Troubleshooting](docs/TROUBLESHOOTING.md) covers host compatibility,
installation, and parameter errors.
[Release validation](docs/RELEASE_VALIDATION.md) records the completed build,
deployment, numerical-example, and downloaded-checkout checks.

## Choose a calculation

| Algorithm | Statistical object | Main controls |
|---|---|---|
| `HOSC_TEMPLATE` | Joint cumulant of values at translated template locations | Integer direction vectors, lag steps, occurrences, Bell/AVL |
| `HOSC_WINDOW` | Local univariate moment or cumulant in each moving window | Order, window size, move step, normalized weights |

The exact template mode supports two- and three-dimensional grids and total
cumulant orders up to eight. It uses unweighted, complete translated templates;
the output axes enumerate lag combinations.

The window mode supports orders one through eight. Available weighting schemes
are `uniform`, `coordinate_gaussian`, and `soft_high_porosity`. Optional outputs
include standardized cumulants, support, effective sample size, validity flags,
and double-precision audit files. The output grid represents window centers.
Bell and AVL select the unweighted evaluators; enabled window weighting uses
the shared weighted recurrence.

## Code, validation, and paper data

- [Repository contents](docs/REPOSITORY_CONTENTS.md): package and source layout.
- [Numerical contract](docs/NUMERICAL_CONTRACT.md) and
  [parameter mapping](docs/PARAMETER_MAPPING.md): statistical definitions and
  input/output rules.
- [Tests](tests/README.md) and [manual checks](docs/MANUAL_TESTING.md): expected
  values and additional integration checks.
- [Standalone implementation](standalone/README.md): C++17 command-line
  calculations without SGeMS, with optional Python experiment workflows.
- [Paper data](Paper-data/README.txt): compact results, figure inputs, provenance,
  and upstream reservoir-data instructions.

The synthetic examples verify software behavior. The scientific results are
provided separately under `Paper-data/current_results/2026-10-12/`. The
[current Supplementary Data S1](Paper-data/supplementary_data_s1/README.md)
adds the manuscript table mappings and matched-window evidence for Fig. 8,
Table 4, Supplementary Table S11, and Fig. S16.

Source rebuilding is optional. Plugin source builds require a matching
VC7.1/SGeMS/Qt3 development environment; the standalone programs use a modern
C++17 toolchain. See [the build guide](docs/BUILD_FROM_SOURCE.md).
Dependency and data notices are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
