# HOSC plugins for 64-bit SGeMS

Template-based joint moments and cumulants, and moving-window local statistics, calculated directly inside **SGeMS x64 Beta**. Both modules accept Cartesian grids and irregular point sets.

The supplied Windows plugins target **x64 / Visual C++ 2010 / Qt 4.7.3**. They include the current Qt4 interfaces, external observation weights, and per-run time and memory records. No separate calculation executable or Python installation is required to use the plugins.

## Install and start

1. Close SGeMS and run `Install.cmd "D:\SGeMS-x64-Beta"`.
2. Start `sgems-x64.exe`. Under **Utilities**, select **HOSC_TEMPLATE** or **HOSC_WINDOW**.
3. Load an input object, then use **Load Parameters** for a matching file in [examples/Parameters](examples/Parameters).
4. Choose a persistent **Record directory** before running if full-precision results and performance records are needed.

The installer checks the four plugin files and the supported host libraries, then backs up any replaced plugin files. It does not install SGeMS or alter its other libraries. See [Quick start](docs/QUICKSTART.md).

## What is included

| Directory | Contents |
|---|---|
| [bin/win64](bin/win64) | Ready-to-install DLL/UI pairs and checksums |
| [src](src) and [ui](ui) | In-process numerical core, x64 host adapters, two plugin entry points and Qt4 forms |
| [examples](examples/README.md) | Two small input objects and six focused, loadable parameter examples |
| [paper/Data](paper/Data) | Sections 4.2–4.5 and supplementary weighted-attribute comparison: 138 parameters, 17 result grids and 14 CSV tables |
| [docs](docs/REPOSITORY_CONTENTS.md) | Usage, build environment, weighting and performance definitions |
| [tests](tests/README.md) | Native SGeMS integration, Qt4 interface and numerical regression tests |
| [tools](tools) | x64 build, installation, example checks and parameter-path preparation |

## Calculations

**HOSC_TEMPLATE** estimates joint moments and cumulants from complete matched data events. Its output axes describe lag configurations. **HOSC_WINDOW** maps local univariate statistics to window centers in the original spatial domain.

Both support uniform weighting, cell declustering, polygonal declustering, global kriging declustering and externally supplied observation weights. See the [interface guide](docs/UI_GUIDE.md), [external weights](docs/EXTERNAL_WEIGHTS.md) and [numerical contract](docs/NUMERICAL_CONTRACT.md).

Every executed run writes `performance.csv` and `performance.json` alongside the full-precision output. Runtime and sampled working-set measurements are defined in [Performance records](docs/PERFORMANCE.md). Native SGeMS display properties use float32; use the full-precision files for numerical comparisons.

## Build and reproduce

- [Build the x64 plugins](docs/BUILD_FROM_SOURCE.md) using a compatible external SDK; [SDK sources](docs/SDK.md).
- [Run the small checks](tests/README.md); [validation scope](docs/VALIDATION.md).
- Check package completeness without running calculations: `python tools/check_repository.py`.
- [Use the current paper data](paper/README.md). The parameter files do not load datasets automatically.
- [Third-party notices](THIRD_PARTY_NOTICES.md).

## Developer reference

The optional [numerical reference](dev/reference/README.md) is isolated under `dev/reference/` for regression comparisons. It is not required to install, build or run the native plugins and is not the source of the current paper performance measurements. Old 32-bit packages, Qt3 scripts, obsolete examples and historical experiment output are outside the active distribution.
