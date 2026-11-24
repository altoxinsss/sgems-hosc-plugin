# Validation scope

The current release is the native x64 DLL/UI pair in `bin/win64`. Production numerical code is in `src/native`; `tests/native_core_sha256.csv` records its imported baseline. The optional `dev/reference` implementation is separate and is not compiled into, or invoked by, the plugin build.

## Checks available in this repository

| Check | Establishes | Does not establish |
|---|---|---|
| `tools/check_repository.py` | Source/release/data hashes, x64 packaging, Qt4 forms, parameter selections and local links | Runtime or numerical correctness |
| `tools/run_examples.py` | Six current examples run inside SGeMS, agree with independent small-sample calculations and save performance records | Full paper-scale performance |
| `tests/current_host_checks.py` | Grid/point geometry, support, Bell/AVL agreement, missing values and small declustering cases | Every possible input or model |
| `qt4_ui_smoke.exe` | Real-host selector, form loading, saved selections, feedback and Load/Run/Save | Every display scaling or screen configuration |
| `run_performance_checks.exe` | Per-run timing and sampled memory recorder behavior | Exact instantaneous peak or comparable memory without controlled host state |
| Optional reference suites | Agreement of DLL outputs with the independent developer implementation on matched inputs | Replacement of the native host performance measurements |

Commands and runtime prerequisites are in [tests](../tests/README.md). Test reports are written to the chosen local output directories. No historical report is presented as a fresh execution.

## Paper evidence

The compact [paper data](../paper/README.md) supplies current inputs, loadable parameters, result grids and full-precision CSV records. Native SGeMS display properties use float32; numerical agreement must be checked using the full-precision records. Supplementary weighted-attribute results reuse the matching fifth-order fields as documented in the paper data guide. Software smoke tests do not replace the archived scientific results or their experimental interpretation.
