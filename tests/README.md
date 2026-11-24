# Tests for the native SGeMS x64 plugins

The primary checks load the DLLs inside a real SGeMS process. They do not replace plugin runs with the optional reference executables. Run commands from the repository root with Python 3.10 or later. Generated files belong in ignored `build/` and `output/` directories, not in the paper evidence package.

## 1. Check the distribution without running SGeMS

```powershell
python tools/check_repository.py
python tests/repository_tools_checks.py
.\tools\install_x64.ps1 -TargetRoot "D:\SGeMS-x64-Beta" -VerifyOnly
```

These verify the numerical-source baseline, release hashes, x64 DLL headers, UI copies, paper-data hashes, evaluator selections and documentation links. Parameter-copy tests also confirm that only paths change and that existing output is not overwritten. They do not demonstrate numerical correctness on their own.

## 2. Run the six introductory examples

Prepare an isolated runtime in a **new directory**. The source installation is read but never modified.

```powershell
python tools/stage_test_host.py --host "D:\SGeMS-x64-Beta" --output "D:\HOSC-tests\host"
python tools/run_examples.py --host "D:\HOSC-tests\host" --out "D:\HOSC-tests\examples"
```

The runner checks all six bundled parameter files against independent fourth-order partition calculations, saved native results and successful 64-bit performance records. The output path must be new to avoid stale results.

## 3. Geometry, support and evaluator regressions

```powershell
python tests/current_host_checks.py --host "D:\HOSC-tests\host" --out "D:\HOSC-tests\core"
```

This suite covers grid and point input, exact and tolerance template searches, raw moments, Bell/AVL selections, missing values, origin/spacing, output slices, support, standardization and the three declustering methods. It includes independent partition calculations. Use a fresh output directory for each run.

## 4. Build and check the Qt4 interface and run recorder

```powershell
python tools/build_x64.py --sdk-root "D:\Development\DLL\DLL\WithSGeMS64\sdk" --with-tests
```

After step 3, run the interface test with the same isolated host and its successful fixtures:

```powershell
python tools/run_ui_checks.py --host "D:\HOSC-tests\host" --suite "D:\HOSC-tests\core" --out "D:\HOSC-tests\ui"
```

This invokes `build/x64/qt4_ui_smoke.exe`, loads the actual host selector and Qt4 forms, checks saved evaluator/weight selections, and runs and saves outputs through the control panel. `--exe` selects another build directory. The runner sets the host DLL and embedded Python paths for the child process only. The form is placed off-screen by the test. Saved results are compared with the independently checked host fixtures.

`build/x64/run_performance_checks.exe <existing-output-directory>` checks per-run timer, sampled memory peak, reset versus lifetime peak, success/failure records and monitor-handle cleanup. It uses a temporary 64 MiB allocation, not a paper benchmark.

To test newly built DLLs, stage another fresh host with `--dll-dir build/x64`. Builds do not overwrite the distributed or installed plugins.

## Optional developer-reference comparisons

The following tests also need the executables built from [dev/reference](../dev/reference/README.md):

- `external_weight_host_checks.py`: external weights, all-one/scaled weights, missing data and invalid weights.
- `window_declustering_host_checks.py`: window weighting and geometry against the independent reference.
- `window_declustering_edge_checks.py`: edge cases for the same comparison.
- `native_reference_checks.py`: comparison helpers and extended historical fixture replay; its `--reference` argument requires an earlier run archive, which is not included in the minimal distribution.

Use each script's `--help` for its required paths. These are development checks, not additional software required by SGeMS users. See [validation scope](../docs/VALIDATION.md).
