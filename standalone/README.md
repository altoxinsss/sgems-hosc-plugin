# Standalone spatial cumulants

The C++17 command-line programs calculate template-based spatial cumulants
(`template_cli`) and moving-window marginal cumulants (`movewindow_cli`) without
installing SGeMS. The numerical C++ sources are copied from the maintained
`WithoutSGeMS` implementation. Bell, AVL and weighted-moment implementations,
compatibility settings, native checks, and small synthetic grid fixtures are
included.

## Build and check

The verified toolchain is Windows with MSVC and CMake 3.20 or newer. A C++17
compiler is required. The retained legacy cores use the MSVC `_isnan` function;
other compiler environments need a corresponding compatibility definition and
have not been validated for this publication snapshot. Run from this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

Visual Studio generators place the executables in `build/Release`; single
configuration generators normally place them in `build`. CTest generates its
input files automatically and checks weighted cumulants, the moving-window AVL
core, both CLIs, and signed integer template geometry.

A small moving-window example is available immediately after building:

```powershell
build/Release/movewindow_cli.exe build/tests/movewindow_cli_smoke.ini
```

With a single configuration generator, use `build/movewindow_cli` instead.
Results are written under `build/test_output`. Input grid and output paths in
an INI file are resolved relative to that file. See
[SGeMS compatibility](SGEMS_COMPATIBILITY.md),
[weighting schemes](WEIGHTING_SCHEMES.md), and
[configuration fragments](configs/weighting/) for the available options.

## Experiment code

The publication subset contains the direct-statistics Sections 4.4–4.6 workflow
dated 2026-10-11, its 2026-10-12 conditional-window supplement, and the order 4–8
scalability workflow. It includes preprocessing,
directional surrogate analysis, controlled-distribution checks, weighted local
maps, spatial block bootstrap intervals, figures, reports, and evidence checksums.
The `r1_5_incremental_value_20260924` directory name is retained because the active
workflow imports shared numerical and reporting modules from that location.

Python 3.11 or newer is required:

```sh
python -m pip install -r requirements.txt
python -m pytest experiments -q
```

Tests use `build/Release` on Windows and `build` elsewhere. For another generator
or location, set `WITHOUT_SGEMS_BUILD_DIR` to the directory containing the built
executables. This setting is also supplied automatically by the direct runner.

See [the reproduction guide](experiments/README.md) for source data, local path
configuration, and execution commands. The native and Python checks use small
fixtures for code acceptance. Formal scientific results require a completed
formal run and its sealed evidence bundle.

Generated builds, local configurations, external raw data, and experiment outputs
belong in ignored directories. The publication includes source code and compact
fixtures; generated large maps and historical run directories are not required
to build or test the programs.
